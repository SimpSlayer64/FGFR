"""Five-seed summary (seeds 11, 22, 33, 44, 55; exhaustiveness 8; box 30) and the direction call from the 95 percent t-interval.
Reads run metadata and rmsd.csv; writes replicate_scores.csv (one row per run) and replicate_summary.csv (one row per compound). No docking."""
import os, sys, csv, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_lib as L
log = L.Log('p2_replicates')
comps = list(csv.DictReader(open(f'{L.REV}/compounds.csv'))); order = [c['compound'] for c in comps]
paper = {r['compound']: r for r in csv.DictReader(open(f'{L.REV}/paper_values.csv'))}
V2 = {'receptor_fgfr3_v2.pdbqt': 'fgfr3', 'receptor_fgfr1_v2.pdbqt': 'fgfr1'}
S = {}
for m in L.all_meta():
    if m['group'] != 'replicates' or int(m['seed']) not in L.SEEDS5 or m['receptor'] not in V2 or int(m['box']) != 30 or int(m['exh']) != 8: continue
    c = os.path.basename(m['ligand'])[:-len('_rebuilt.pdbqt')]
    S[(c, V2[m['receptor']], int(m['seed']))] = m
# core-RMSD (copy A, model 1) per run
R = {}
if os.path.exists(f'{L.REV}/rmsd.csv'):
    for r in csv.DictReader(open(f'{L.REV}/rmsd.csv')):
        if r['stage'] == '2b_replicates' and r['model'] == '1' and r['crystal_copy'] == 'A' and r['status'] == 'ok':
            R[(r['compound'], r['receptor'], int(r['seed']))] = float(r['rmsd'])
runs = []; summ = []
for c in order:
    ok = all((c, t, s) in S for t in ('fgfr3', 'fgfr1') for s in L.SEEDS5)
    if not ok: log.p(f'  incomplete: {c}'); continue
    si = []; f3s = []; f1s = []
    for s in L.SEEDS5:
        a, b = float(S[(c, 'fgfr3', s)]['best_score']), float(S[(c, 'fgfr1', s)]['best_score']); si.append(a - b); f3s.append(a); f1s.append(b)
        runs.append(dict(compound=c, in_paper_20='yes' if c in paper else 'NO (hybrid1)', seed=s, fgfr3_score=a, fgfr1_score=b, si=round(a - b, 3),
                         core_rmsd_fgfr3=R.get((c, 'fgfr3', s), ''), core_rmsd_fgfr1=R.get((c, 'fgfr1', s), ''),
                         fgfr3_job=S[(c, 'fgfr3', s)]['id'], fgfr1_job=S[(c, 'fgfr1', s)]['id'], vina_version=S[(c, 'fgfr3', s)]['vina_version'], timestamp=S[(c, 'fgfr3', s)]['timestamp']))
    m, sd, lo, hi, out = L.rule_outcome(si)
    m3, sd3 = L.mean_sd(f3s); m1, sd1 = L.mean_sd(f1s)
    pv = paper.get(c)
    n_bad3 = sum(1 for s in L.SEEDS5 if R.get((c, 'fgfr3', s), 0) > 3.0); n_bad1 = sum(1 for s in L.SEEDS5 if R.get((c, 'fgfr1', s), 0) > 3.0)
    n_bad3b = sum(1 for s in L.SEEDS5 if R.get((c, 'fgfr3', s), 0) > 2.0); n_bad1b = sum(1 for s in L.SEEDS5 if R.get((c, 'fgfr1', s), 0) > 2.0)
    summ.append(dict(compound=c, in_paper_20='yes' if pv else 'NO (hybrid1)', mean_fgfr3=round(m3, 3), sd_fgfr3=round(sd3, 3), mean_fgfr1=round(m1, 3), sd_fgfr1=round(sd1, 3),
                     si_per_seed=' '.join(f'{x:+.2f}' for x in si), mean_si=round(m, 3), sd_si=round(sd, 3), ci_lo=round(lo, 3), ci_hi=round(hi, 3), rule_outcome=out,
                     margin=round(abs(m) - L.T975_DF4 * sd / math.sqrt(5), 3),
                     paper_fgfr3=pv['fgfr3_corr'] if pv else '', paper_fgfr1=pv['fgfr1_corr'] if pv else '', paper_si=pv['si_corr'] if pv else '',
                     paper_sign_matches_outcome=('' if not pv else ('yes' if (float(pv['si_corr']) < 0 and out == 'FGFR3-preferring') or (float(pv['si_corr']) > 0 and out == 'FGFR1-preferring') else ('indistinguishable' if out == 'indistinguishable' else 'NO'))),
                     seeds_core_rmsd_gt_2_fgfr3=n_bad3b, seeds_core_rmsd_gt_3_fgfr3=n_bad3, seeds_core_rmsd_gt_2_fgfr1=n_bad1b, seeds_core_rmsd_gt_3_fgfr1=n_bad1))
for fn, rows in (('replicate_scores.csv', runs), ('replicate_summary.csv', summ)):
    if rows:
        with open(f'{L.REV}/{fn}', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
log.p(f'{len(summ)} compounds complete, {len(runs)} runs')
log.p(f'{"compound":16}{"F3":>8}{"F1":>8}{"meanSI":>8}{"SD":>7}{"95% interval":>18}  outcome (paper SI)')
for r in sorted(summ, key=lambda r: r['mean_si']):
    log.p(f"{r['compound']:16}{r['mean_fgfr3']:8.3f}{r['mean_fgfr1']:8.3f}{r['mean_si']:+8.3f}{r['sd_si']:7.3f}  [{r['ci_lo']:+.3f},{r['ci_hi']:+.3f}]  {r['rule_outcome']:18} paper {r['paper_si']}")
cnt = {}
for r in summ: cnt[r['rule_outcome']] = cnt.get(r['rule_outcome'], 0) + 1
log.p('rule outcomes:', cnt)
d = [r for r in summ if r['compound'] == 'dabogratinib']
if d: log.p(f"dabogratinib SI over five seeds: mean {d[0]['mean_si']:+.3f}, SD {d[0]['sd_si']:.3f}")
