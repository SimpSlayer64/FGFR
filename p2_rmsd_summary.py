"""Summaries of rmsd.csv: redock of dabogratinib (2f) and pose counts at 2.0 and 3.0 A (2b). Usage: python p2_rmsd_summary.py"""
import os, sys, csv, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_lib as L
log = L.Log('p2_rmsd_summary')
rows = list(csv.DictReader(open(f'{L.REV}/rmsd.csv')))
def f(x): return float(x) if x not in ('', None) else None
copyA = {'fgfr3': 'A', 'fgfr1': 'A'}
out = []
log.p('== 2f: validation redock, dabogratinib vs crystal ligand (full heavy-atom, symmetry-aware, in place) ==')
for st in ('2f_seed42', '2f_fiveseed'):
    for t in ('fgfr3', 'fgfr1'):
        sub = [r for r in rows if r['stage'] == st and r['receptor'] == t and r['status'] == 'ok']
        if not sub: continue
        for cp in sorted({r['crystal_copy'] for r in sub}):
            m1 = [(r['seed'], f(r['rmsd']), f(r['centroid_dist_heavy'])) for r in sub if r['crystal_copy'] == cp and r['model'] == '1']
            best5 = {}
            for r in sub:
                if r['crystal_copy'] == cp: best5[r['seed']] = min(best5.get(r['seed'], 99), f(r['rmsd']))
            log.p(f'  {st} {t} copy {cp}: model-1 RMSD ' + ', '.join(f'seed{s}={v:.2f}' for s, v, _ in m1) + ' | centroid ' + ', '.join(f'{c:.2f}' for _, _, c in m1) + ' | best-of-5-models ' + ', '.join(f'{v:.2f}' for v in best5.values()))
            out.append(dict(block='2f', stage=st, receptor=t, copy=cp, model1_rmsd=';'.join(f'{v:.3f}' for _, v, _ in m1), centroid=';'.join(f'{c:.3f}' for _, _, c in m1), best5=';'.join(f'{v:.3f}' for v in best5.values())))
log.p('\n== 2b: counts of model-1 poses by RMSD (core for analogs, same rule both receptors) ==')
def counts(sub, label, t):
    vals = [v for v in sub if v is not None]; n = len(sub)
    d = dict(block='2b', set=label, receptor=t, n_poses=n, n_with_rmsd=len(vals), le_2p0=sum(v <= 2.0 for v in vals), le_3p0=sum(v <= 3.0 for v in vals), gt_3p0=sum(v > 3.0 for v in vals), core_too_small=n - len(vals))
    log.p(f"  {label:34}{t:6} n={n:3d} with RMSD {len(vals):3d}: <=2.0 A {d['le_2p0']:3d} | <=3.0 A {d['le_3p0']:3d} | >3.0 A {d['gt_3p0']:3d} | core_too_small {d['core_too_small']}")
    out.append(d)
for st, lab in (('2b_seed42', 'seed-42 set (paper poses)'), ('2b_replicates', 'five-seed pooled set')):
    for t in ('fgfr3', 'fgfr1'):
        sub = [r for r in rows if r['stage'] == st and r['receptor'] == t and r['model'] == '1']
        if not sub: continue
        byp = {}
        for r in sub: byp.setdefault((r['compound'], r['seed']), {})[r['crystal_copy']] = r
        for scope in ('copy A', 'nearest copy'):
            vals = []
            for k, d in byp.items():
                v = [f(x['rmsd']) for x in d.values() if x['status'] == 'ok' and f(x['rmsd']) is not None]
                if scope == 'copy A': vals.append(f(d['A']['rmsd']) if 'A' in d and d['A']['status'] == 'ok' else None)
                else: vals.append(min(v) if v else None)
            counts(vals, f'{lab}, {scope}', t)
        if st == '2b_replicates':
            comp = {}
            for (c, s), d in byp.items():
                if d['A']['status'] == 'ok': comp.setdefault(c, []).append(f(d['A']['rmsd']))
            med = [statistics.median(v) for v in comp.values()]
            counts(med + [None] * (len({k[0] for k in byp}) - len(med)), 'compound median over seeds, copy A', t)
log.p('\n== cross-tab: seed-42 core RMSD (copy A, model 1) vs the paper centroid-drift flag (>= 3.3 A, pose_drift_corr.csv) ==')
drift = {r['compound']: r for r in csv.DictReader(open(f'{L.REV}/pose_drift_corr.csv'))}
for t, col in (('fgfr3', 'drift3_allatom'), ('fgfr1', 'drift1_allatom')):
    tab = {}
    for r in rows:
        if r['stage'] == '2b_seed42' and r['receptor'] == t and r['model'] == '1' and r['crystal_copy'] == 'A' and r['status'] == 'ok':
            cflag = float(drift[r['compound']][col]) >= 3.3
            rf = 'rmsd<=2' if f(r['rmsd']) <= 2.0 else ('rmsd<=3' if f(r['rmsd']) <= 3.0 else 'rmsd>3')
            tab.setdefault((cflag, rf), []).append(r['compound'])
    for k in sorted(tab): log.p(f'  {t}: centroid flag (>=3.3 A) {str(k[0]):5} & {k[1]:8}: {len(tab[k])}  {tab[k] if len(tab[k]) <= 8 else ""}')
with open(f'{L.REV}/rmsd_summary.csv', 'w', newline='') as fh:
    keys = sorted({k for d in out for k in d}); w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(out)
