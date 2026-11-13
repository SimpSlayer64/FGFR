"""Ensemble analysis: fraction of FGFR3-FGFR1 receptor pairs with a negative SI, distribution, rank of the reference pair, 4 A site filter.
Reads runs/ensemble/*.json and the pose files; writes phase3b/ensemble/*.csv and ENSEMBLE_RESULT.txt.
Usage: PYTHONDONTWRITEBYTECODE=1 python p3b_analyze.py"""
import sys, os, csv, json, glob, statistics as stt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_lib as L
import numpy as np
B = f'{L.REV}/phase3b'; E = f'{B}/ensemble'
log = L.Log('p3b_analyze')
CUT = 4.0
R = {r['entry']: r for r in csv.DictReader(open(f'{E}/receptors.csv'))}
RL = {r['entry']: r for r in csv.DictReader(open(f'{B}/receptor_list.csv'))}
def centroid(pose):
    xyz = []
    for l in open(pose):
        if l.startswith('ENDMDL'): break
        if l.startswith(('ATOM', 'HETATM')) and l[77:79].strip() not in ('HD', 'H'): xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(xyz).mean(0), len(xyz)
runs = []
for p in sorted(glob.glob(f'{L.RUNS}/ensemble/*.json')):
    m = json.load(open(p)); entry = m['tag'].split()[-1]; r = R[entry]
    c, n = centroid(f"{L.REV}/{m['pose_file']}"); ref = np.array([float(x) for x in r['ref9CD7_ligand_centroid_in_receptor_frame'].split()])
    runs.append(dict(entry=entry, kind=r['kind'], isoform=r['isoform'], seed=int(m['seed']), score=float(m['best_score']), top_pose_heavy_atoms=n, centroid_distance_A=round(float(np.linalg.norm(c - ref)), 3), seconds=m['seconds'], job_id=m['id']))
with open(f'{E}/ensemble_per_run.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(runs[0].keys())); w.writeheader(); w.writerows(runs)
ens_ids = [e for e, r in R.items() if r['kind'] == 'ensemble']; need = {(e, s) for e in R for s in (11, 22, 33)}
have = {(r['entry'], r['seed']) for r in runs}; log.p(f'runs found {len(runs)} of {len(need)} required; missing {sorted(need - have)}')
assert not (need - have), 'incomplete'
def receptor_scores(filt, exclude=()):
    out = {}
    for e in R:
        rr = [r for r in runs if r['entry'] == e and (not filt or r['centroid_distance_A'] <= CUT)]
        if rr: out[e] = dict(score=float(np.mean([r['score'] for r in rr])), n_runs=len(rr), isoform=R[e]['isoform'], kind=R[e]['kind'])
    return out
def analyse(name, filt, exclude=()):
    rs = receptor_scores(filt); ens = {e: v for e, v in rs.items() if v['kind'] == 'ensemble' and e not in exclude}
    f3 = {e: v['score'] for e, v in ens.items() if v['isoform'] == 'FGFR3'}; f1 = {e: v['score'] for e, v in ens.items() if v['isoform'] == 'FGFR1'}
    pairs = [(a, b, f3[a] - f1[b]) for a in f3 for b in f1]; si = np.array([p[2] for p in pairs])
    res = dict(variant=name, n_fgfr3=len(f3), n_fgfr1=len(f1), n_pairs=len(si), fraction_negative=float((si < 0).mean()), median=float(np.median(si)), p5=float(np.percentile(si, 5)), p95=float(np.percentile(si, 95)), mean=float(si.mean()), min=float(si.min()), max=float(si.max()))
    for iso, d in (('FGFR3', f3), ('FGFR1', f1)):
        v = np.array(list(d.values())); res[f'{iso}_n'] = len(v); res[f'{iso}_min'] = float(v.min()); res[f'{iso}_p5'] = float(np.percentile(v, 5)); res[f'{iso}_median'] = float(np.median(v)); res[f'{iso}_mean'] = float(v.mean()); res[f'{iso}_p95'] = float(np.percentile(v, 95)); res[f'{iso}_max'] = float(v.max()); res[f'{iso}_sd'] = float(v.std(ddof=1))
    # reference pair
    if '9CD7' in rs and '9CD5' in rs:
        rsi = rs['9CD7']['score'] - rs['9CD5']['score']; below = int((si < rsi).sum()); res.update(ref_si=float(rsi), ref_fgfr3=rs['9CD7']['score'], ref_fgfr1=rs['9CD5']['score'], ref_n_runs=f"{rs['9CD7']['n_runs']}/{rs['9CD5']['n_runs']}", ref_rank_among_N_plus_1=below + 1, ref_N_plus_1=len(si) + 1, ref_fraction_of_pairs_below=below / len(si))
    # per receptor fractions
    per3 = {a: float(np.mean([p[2] < 0 for p in pairs if p[0] == a])) for a in f3}; per1 = {b: float(np.mean([p[2] < 0 for p in pairs if p[1] == b])) for b in f1}
    rng = np.random.default_rng(20261011); k3, k1 = list(f3), list(f1); M = np.array([[f3[a] - f1[b] for b in k1] for a in k3]); fr = []
    for _ in range(10000):
        i = rng.integers(0, len(k3), len(k3)); j = rng.integers(0, len(k1), len(k1)); fr.append((M[np.ix_(i, j)] < 0).mean())
    res['bootstrap_fraction_negative_2.5'] = float(np.percentile(fr, 2.5)); res['bootstrap_fraction_negative_97.5'] = float(np.percentile(fr, 97.5))
    F = res['fraction_negative']; res['reading'] = 'sign belongs to the setup (F >= 0.80)' if F >= 0.8 else ('the 9CD7/9CD5 pair is special (F <= 0.20)' if F <= 0.2 else 'no reading (0.20 < F < 0.80); observed fraction reported')
    return res, pairs, per3, per1, ens
variants = [('A_all_runs', False, ()), ('B_without_runs_gt_4A', True, ()), ('S1_A_without_pocket_difference_receptors', False, ('8UDV', '5VND', '6P69')), ('S1_B_without_pocket_difference_receptors', True, ('8UDV', '5VND', '6P69'))]
allres = []; detail = {}
for name, filt, ex in variants:
    res, pairs, per3, per1, ens = analyse(name, filt, ex); allres.append(res); detail[name] = (pairs, per3, per1, ens)
keys = list(dict.fromkeys(k for r in allres for k in r))
with open(f'{E}/ensemble_summary.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(allres)
with open(f'{E}/ensemble_pairs.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['variant', 'fgfr3_entry', 'fgfr1_entry', 'SI']); 
    for name in detail:
        for a, b, s in detail[name][0]: w.writerow([name, a, b, round(s, 4)])
with open(f'{E}/ensemble_receptor_scores.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['variant', 'entry', 'isoform', 'n_runs_used', 'mean_score', 'fraction_of_pairs_negative'])
    for name in detail:
        pairs, per3, per1, ens = detail[name]
        for e, v in ens.items(): w.writerow([name, e, v['isoform'], v['n_runs'], round(v['score'], 3), round((per3 if v['isoform'] == 'FGFR3' else per1)[e], 3)])
# ---------- text report
t = []
P = t.append
P('ENSEMBLE RESULT (computed by p3b_analyze.py)')
dist = np.array([r['centroid_distance_A'] for r in runs if r['kind'] == 'ensemble'])
P(f"top-pose centroid distance to the superposed 9CD7 ligand centroid, ensemble runs (n={len(dist)}): min {dist.min():.2f} median {np.median(dist):.2f} max {dist.max():.2f} A; runs > 4 A: {int((dist > CUT).sum())} ({(dist > CUT).mean():.1%}); by isoform > 4 A: FGFR1 {sum(1 for r in runs if r['kind']=='ensemble' and r['isoform']=='FGFR1' and r['centroid_distance_A']>CUT)}/{sum(1 for r in runs if r['kind']=='ensemble' and r['isoform']=='FGFR1')}, FGFR3 {sum(1 for r in runs if r['kind']=='ensemble' and r['isoform']=='FGFR3' and r['centroid_distance_A']>CUT)}/{sum(1 for r in runs if r['kind']=='ensemble' and r['isoform']=='FGFR3')}")
P('reference runs (9CD7, 9CD5): ' + '; '.join(f"{r['entry']} s{r['seed']} score {r['score']} dist {r['centroid_distance_A']} A" for r in runs if r['kind'] == 'reference'))
for r in allres:
    P(f"\n[{r['variant']}] FGFR3 receptors {r['n_fgfr3']}, FGFR1 receptors {r['n_fgfr1']}, pairs {r['n_pairs']}")
    P(f"  fraction negative {r['fraction_negative']:.3f} | SI median {r['median']:+.3f}, 5th {r['p5']:+.3f}, 95th {r['p95']:+.3f}, mean {r['mean']:+.3f}, range [{r['min']:+.3f}, {r['max']:+.3f}] | receptor-level bootstrap of the fraction 95% [{r['bootstrap_fraction_negative_2.5']:.3f}, {r['bootstrap_fraction_negative_97.5']:.3f}]")
    for iso in ('FGFR3', 'FGFR1'): P(f"  {iso} receptor scores (n={r[iso+'_n']}): min {r[iso+'_min']:.3f} 5th {r[iso+'_p5']:.3f} median {r[iso+'_median']:.3f} mean {r[iso+'_mean']:.3f} 95th {r[iso+'_p95']:.3f} max {r[iso+'_max']:.3f} SD {r[iso+'_sd']:.3f}")
    P(f"  9CD7/9CD5 pair: FGFR3 {r['ref_fgfr3']:.3f} FGFR1 {r['ref_fgfr1']:.3f} SI {r['ref_si']:+.3f} (runs used {r['ref_n_runs']}); rank {r['ref_rank_among_N_plus_1']} of {r['ref_N_plus_1']} (1 = most negative); fraction of ensemble pairs below it {r['ref_fraction_of_pairs_below']:.3f}")
    P(f"  READING: {r['reading']}")
    pairs, per3, per1, ens = detail[r['variant']]
    P('  per FGFR3 receptor: fraction of FGFR1 partners with negative SI and score: ' + '; '.join(f"{e} {per3[e]:.2f} ({ens[e]['score']:.3f})" for e in per3))
    P('  per FGFR1 receptor: fraction of FGFR3 partners with negative SI: ' + ' '.join(f"{e}:{per1[e]:.1f}" for e in sorted(per1)))
open(f'{E}/ENSEMBLE_RESULT.txt', 'w').write('\n'.join(t) + '\n'); log.p('\n'.join(t))
