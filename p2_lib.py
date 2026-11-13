"""Shared helpers: paths, logging, docking job runner (Vina Python API), manifests and summary statistics.
Paths come from the environment: FGFR3_DOCK (folder with the structures and prepared receptors, default: current folder)
and FGFR3_REV (working folder for runs, logs and manifests, default: $FGFR3_DOCK/work)."""
import os, sys, json, csv, time, glob, platform, datetime, hashlib, math

DOCK = os.environ.get('FGFR3_DOCK', os.getcwd())          # folder holding 9CD7.cif, 9CD5.cif and the prepared receptors
REV = os.environ.get('FGFR3_REV', f'{DOCK}/work')    # working folder for runs, logs, manifests and ensemble outputs
RUNS = f'{REV}/runs'
LOGS = f'{REV}/logs'
T975_DF4 = 2.776     # Student t, 95 percent two-sided, 4 degrees of freedom
SEEDS5 = [11, 22, 33, 44, 55]

def versions_header(script, extra=''):
    import rdkit, numpy, meeko, gemmi
    from vina import __version__ as vv
    return (f'# {script}  {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}  args={sys.argv[1:]} {extra}\n'
            f'# python {platform.python_version()} | vina {vv} | meeko {meeko.__version__} | rdkit {rdkit.__version__} | gemmi {gemmi.__version__} | numpy {numpy.__version__}')

class Log:
    def __init__(self, name, extra=''):
        os.makedirs(LOGS, exist_ok=True)
        self.f = open(f'{LOGS}/{name}.log', 'a')
        self.p(versions_header(name, extra))
    def p(self, *a):
        s = ' '.join(str(x) for x in a); print(s, flush=True); self.f.write(s + '\n'); self.f.flush()

def sha256(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()

def box_centres():
    j = json.load(open(f'{DOCK}/box_centers.json'))
    return {'fgfr3': j['fgfr3'], 'fgfr1': j['fgfr1']}

# ---------------------------------------------------------------- docking jobs
JOB_FIELDS = ['id', 'group', 'ligand', 'receptor', 'cx', 'cy', 'cz', 'box', 'exh', 'seed', 'cpu', 'n_poses', 'tag']
def job_id(ligand, receptor, centre, box, exh, seed, conf):
    return (f'{os.path.basename(ligand)[:-6]}__{os.path.basename(receptor)[:-6]}__c{centre[0]:.2f}_{centre[1]:.2f}_{centre[2]:.2f}'
            f'__b{box}__e{exh}__s{seed}__{conf}')
def make_job(group, ligand, receptor, centre, box=30, exh=8, seed=42, cpu=0, conf='conf0', tag=''):
    """ligand/receptor are paths relative to DOCK or absolute."""
    return dict(id=job_id(ligand, receptor, centre, box, exh, seed, conf), group=group, ligand=ligand, receptor=receptor,
                cx=centre[0], cy=centre[1], cz=centre[2], box=box, exh=exh, seed=seed, cpu=cpu, n_poses=5, tag=tag)
def pose_path(job): return f"{RUNS}/{job['group']}/{job['id']}.pdbqt"
def meta_path(job): return f"{RUNS}/{job['group']}/{job['id']}.json"
def abspath(p): return p if os.path.isabs(p) else f'{DOCK}/{p}'

def run_job(job):
    """Dock one job. Skips (returns 'skip') if both pose and meta already exist. Writes pose to temp then renames."""
    out, meta = pose_path(job), meta_path(job)
    if os.path.exists(out) and os.path.exists(meta):
        return ('skip', job['id'], None)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    from vina import Vina, __version__ as vv
    cpu = int(job['cpu'])
    v = Vina(sf_name='vina', cpu=cpu, seed=int(job['seed']), verbosity=0)
    t0 = time.time()
    v.set_receptor(abspath(job['receptor'])); v.set_ligand_from_file(abspath(job['ligand']))
    v.compute_vina_maps(center=[float(job['cx']), float(job['cy']), float(job['cz'])], box_size=[float(job['box'])] * 3)
    v.dock(exhaustiveness=int(job['exh']), n_poses=int(job['n_poses']))
    tmp = out + '.tmp'
    v.write_poses(tmp, n_poses=int(job['n_poses']), overwrite=True)
    en = v.energies(n_poses=int(job['n_poses']))
    dt = time.time() - t0
    os.replace(tmp, out)
    row = dict(job, vina_version=vv, best_score=round(float(en[0][0]), 3), all_scores=' '.join(f'{e[0]:.3f}' for e in en),
               seconds=round(dt, 2), timestamp=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), pose_file=os.path.relpath(out, REV),
               ligand_md5=hashlib.md5(open(abspath(job['ligand']), 'rb').read()).hexdigest(),
               receptor_md5=hashlib.md5(open(abspath(job['receptor']), 'rb').read()).hexdigest())
    json.dump(row, open(meta + '.tmp', 'w')); os.replace(meta + '.tmp', meta)
    return ('done', job['id'], row)

def run_job_safe(job):
    try: return run_job(job)
    except Exception as e: return ('error', job['id'], f'{type(e).__name__}: {e}')

def read_manifest(path):
    return list(csv.DictReader(open(path)))
def write_manifest(path, jobs):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=JOB_FIELDS); w.writeheader()
        for j in jobs: w.writerow({k: j[k] for k in JOB_FIELDS})

def all_meta():
    rows = []
    for p in sorted(glob.glob(f'{RUNS}/*/*.json')):
        try: rows.append(json.load(open(p)))
        except Exception: pass
    return rows
def rebuild_docking_runs_csv():
    rows = all_meta()
    if not rows: return 0
    keys = ['id', 'group', 'ligand', 'ligand_md5', 'receptor', 'receptor_md5', 'cx', 'cy', 'cz', 'box', 'exh', 'seed', 'cpu', 'n_poses', 'vina_version',
            'best_score', 'all_scores', 'seconds', 'timestamp', 'pose_file', 'tag']
    with open(f'{REV}/docking_runs.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore'); w.writeheader(); w.writerows(rows)
    return len(rows)

# ---------------------------------------------------------------- statistics
def mean_sd(x):
    n = len(x); m = sum(x) / n
    sd = math.sqrt(sum((a - m) ** 2 for a in x) / (n - 1)) if n > 1 else float('nan')
    return m, sd
def rule_outcome(si_list):
    """Direction call: 95 percent t-interval (t = 2.776) of the mean SI over five seeds; SD (n-1) of that compound's own runs."""
    m, sd = mean_sd(si_list); hw = T975_DF4 * sd / math.sqrt(len(si_list))
    lo, hi = m - hw, m + hw
    out = 'FGFR3-preferring' if hi < 0 else ('FGFR1-preferring' if lo > 0 else 'indistinguishable')
    return m, sd, lo, hi, out
