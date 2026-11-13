"""Docks dabogratinib into the reference pair (9CD7, 9CD5) and the receptor ensemble. Resumable; one docking at a time; cpu=0.
The reference runs are docked first and must reproduce the stored replicate scores to 0.005. Usage: PYTHONDONTWRITEBYTECODE=1 python p3b_dock.py"""
import sys, os, csv, json, time, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_lib as L
B = f'{L.REV}/phase3b'; E = f'{B}/ensemble'
log = L.Log('p3b_dock')
LIG = 'work/ligands/dabogratinib_rebuilt.pdbqt'; SEEDS = [11, 22, 33]
R = list(csv.DictReader(open(f'{E}/receptors.csv')))
def jobs_for(rows): return [L.make_job('ensemble', LIG, r['receptor_file'], (float(r['cx']), float(r['cy']), float(r['cz'])), box=30, exh=8, seed=s, cpu=0, tag=f"3b {r['kind']} {r['entry']}") for r in rows for s in SEEDS]
ref = jobs_for([r for r in R if r['kind'] == 'reference']); ens = jobs_for([r for r in R if r['kind'] == 'ensemble'])
L.write_manifest(f'{E}/manifest.csv', ref + ens); log.p(f'manifest: {len(ref)} reference + {len(ens)} ensemble = {len(ref) + len(ens)} jobs')
def run(jobs, label):
    t0 = time.time(); done = 0; nerr = 0
    for k, j in enumerate(jobs, 1):
        st, jid, row = L.run_job_safe(j)
        if st == 'error':
            log.p('ERROR (retry once)', jid, row); st, jid, row = L.run_job_safe(j)
            if st == 'error': log.p('ERROR again', jid, row); nerr += 1; continue
        if st == 'done': done += 1
        if k % 10 == 0 or k == len(jobs): log.p(f'  {label}: {k}/{len(jobs)} ({done} newly run) {time.time() - t0:.0f} s elapsed')
        if label == 'ensemble' and done == 20 and st == 'done':
            per = (time.time() - t0) / done; proj = per * (len(jobs)) / 3600; log.p(f'  measured {per:.1f} s per docking over the first 20 runs; projected ensemble {proj:.2f} h (stop if total > 3 h)')
            if proj > 3: log.p('S1: projected runtime exceeds 3 h; stopping'); sys.exit(10)
    return nerr
nerr = run(ref, 'reference')
bad = []
for j in ref:
    new = json.load(open(L.meta_path(j)))['best_score']; old_path = f"{L.RUNS}/replicates/{j['id']}.json"
    old = json.load(open(old_path))['best_score'] if os.path.exists(old_path) else None
    ok = old is not None and abs(new - old) <= 0.005; log.p(f"S2 check {j['id'].split('__')[1]} seed {j['seed']}: new {new} stored replicate {old} {'OK' if ok else 'MISMATCH'}"); 
    if not ok: bad.append(j['id'])
if bad: log.p('S2: reference scores do not reproduce stored replicate; stopping', bad); sys.exit(2)
nerr += run(ens, 'ensemble')
rows = [json.load(open(p)) for p in sorted(glob.glob(f'{L.RUNS}/ensemble/*.json'))]
keys = ['id', 'group', 'ligand', 'ligand_md5', 'receptor', 'receptor_md5', 'cx', 'cy', 'cz', 'box', 'exh', 'seed', 'cpu', 'n_poses', 'vina_version', 'best_score', 'all_scores', 'seconds', 'timestamp', 'pose_file', 'tag']
with open(f'{E}/ensemble_runs.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore'); w.writeheader(); w.writerows(rows)
log.p(f'finished: {len(rows)} run records in ensemble_runs.csv, unresolved errors {nerr}')
sys.exit(4 if nerr else 0)
