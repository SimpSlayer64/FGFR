"""Runs p2_rmsd.py for a stage and merges the result into rmsd.csv. Resumable (a part file that exists is not recomputed).
Usage: PYTHONDONTWRITEBYTECODE=1 python p2_rmsd_stage.py 2f | 2b_seed42 | 2b_replicates
2f: dabogratinib poses (full-molecule RMSD): the seed-42 redock and the five-seed runs (seeds 11-55).
2b_seed42: all 21 seed-42 best-pose sets, shared-core RMSD against the dabogratinib crystal pose (--core-ref).
2b_replicates: every pose file of group 'replicates' (21 compounds x 2 receptors x 5 seeds), shared-core RMSD."""
import os, sys, csv, glob, json, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_lib as L
stage = sys.argv[1]
log = L.Log('p2_rmsd_stage', f'stage={stage}')
DAB = 'C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CC5(C4)CN(S(C)(=O)=O)C5)nc3)c2c1)c1c(Cl)cncc1Cl'
REC = {'fgfr3': ('9CD7.cif', 'A', 'receptor_fgfr3_v2'), 'fgfr1': ('9CD5.cif', 'A', 'receptor_fgfr1_v2')}
names = [r['compound'] for r in csv.DictReader(open(f'{L.REV}/compounds.csv'))]
os.makedirs(f'{L.REV}/rmsd_parts', exist_ok=True)
def replicate_poses(t, compounds=None):
    out = []
    for m in L.all_meta():
        if m['group'] == 'replicates' and m['receptor'] == REC[t][2] + '.pdbqt' and int(m['seed']) in L.SEEDS5:
            c = os.path.basename(m['ligand'])[:-len('_rebuilt.pdbqt')]
            if compounds is None or c in compounds: out.append((c, int(m['seed']), f"{L.REV}/{m['pose_file']}"))
    return sorted(out)
jobs = []   # (receptor, label, poses list [(compound, seed, path)], core_ref)
for t in ('fgfr3', 'fgfr1'):
    if stage == '2f':
        jobs.append((t, '2f_seed42', [('dabogratinib', 42, f'{L.DOCK}/docked_dabogratinib_{t}_corr.pdbqt')], False))
        jobs.append((t, '2f_fiveseed', replicate_poses(t, {'dabogratinib'}), False))
    elif stage == '2b_seed42':
        jobs.append((t, '2b_seed42', [(c, 42, f'{L.DOCK}/docked_{c}_{t}_corr.pdbqt') for c in names], True))
    elif stage == '2b_replicates':
        jobs.append((t, '2b_replicates', replicate_poses(t), True))
allrows = []
for t, label, poses, core in jobs:
    out = f'{L.REV}/rmsd_parts/{label}_{t}_{len(poses)}.csv'
    if not os.path.exists(out):
        cif, ch, _ = REC[t]
        cmd = [sys.executable, f'{L.REV}/p2_rmsd.py', '--cif', f'{L.DOCK}/{cif}', '--chain', ch, '--ligand-resn', 'A1AV2', '--template-smiles', DAB,
               '--poses'] + [p for _, _, p in poses] + ['--out', out] + (['--core-ref'] if core else [])
        r = subprocess.run(cmd, capture_output=True, text=True, env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1'))
        log.p(f'{label} {t}: {len(poses)} pose files, exit {r.returncode}'); open(f'{L.LOGS}/p2_rmsd_{label}_{t}.log', 'w').write(r.stdout + r.stderr)
        if r.returncode != 0: log.p('STOP S3/S4: p2_rmsd.py failed; see logs/p2_rmsd_%s_%s.log' % (label, t)); sys.exit(3)
    else: log.p(f'{label} {t}: part exists, skipped')
    meta = {p: (c, s) for c, s, p in poses}
    for r in csv.DictReader(open(out)):
        c, s = meta[r['pose_file']]; allrows.append(dict(r, stage=label, receptor=t, compound=c, seed=s, pose_file=os.path.relpath(r['pose_file'], L.DOCK) if r['pose_file'].startswith(L.DOCK) else r['pose_file']))
# merge: keep rows of other stages already in rmsd.csv
path = f'{L.REV}/rmsd.csv'; keys = ['stage', 'receptor', 'compound', 'seed', 'pose_file', 'model', 'crystal_copy', 'status', 'rmsd', 'core_atoms', 'n_maps', 'n_pose_heavy', 'centroid_dist_heavy']
old = [r for r in csv.DictReader(open(path))] if os.path.exists(path) else []
labels = {r['stage'] for r in allrows}
rows = [r for r in old if r['stage'] not in labels] + allrows
with open(path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore'); w.writeheader(); w.writerows(rows)
log.p(f'rmsd.csv now {len(rows)} rows (stage {stage}: {len(allrows)} rows)')
