"""Sensitivity test. Usage: PYTHONDONTWRITEBYTECODE=1 python p2_sensitivity.py select | summary
select : choose boundary compounds from replicate_summary.csv, build extra conformers (ETKDGv3 seeds 101, 202) for all subjects, write manifests/d2.csv
summary: compute sensitivity.csv from the finished sensitivity runs plus the base runs at seeds 11, 22, 33."""
import os, sys, csv, json, math, shutil, subprocess, hashlib, importlib.util
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_lib as L
mode = sys.argv[1]
log = L.Log('p2_sensitivity', f'mode={mode}')
SEEDS3 = [11, 22, 33]; MANDATORY = ['dabogratinib', 'v2_iPr_sulfonyl', 'v4_Nme_indazole', 'c3_v3v2', 'c3_v3v2v1']; ALT = MANDATORY[1:]
V2 = {'fgfr3': 'receptor_fgfr3_v2.pdbqt', 'fgfr1': 'receptor_fgfr1_v2.pdbqt'}
bc = L.box_centres()
if mode == 'select':
    summ = list(csv.DictReader(open(f'{L.REV}/replicate_summary.csv')))
    cand = sorted([(abs(float(r['margin'])), r['compound'], float(r['margin'])) for r in summ if abs(float(r['margin'])) <= 0.10 and r['compound'] not in MANDATORY])
    boundary = [c for _, c, _ in cand[:6]]
    log.p('boundary candidates (|margin| <= 0.10, ranked):', [(c, round(m, 3)) for _, c, m in cand]); log.p('boundary compounds kept (cap 6):', boundary)
    subjects = MANDATORY + boundary
    # extra conformers
    spec = importlib.util.spec_from_file_location('build_ligands', f'{L.DOCK}/build_ligands.py'); BL = importlib.util.module_from_spec(spec); sys.dont_write_bytecode = True; spec.loader.exec_module(BL)
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign
    os.makedirs(f'{L.REV}/ligands/conformers', exist_ok=True); mk = shutil.which('mk_prepare_ligand.py')
    for c in subjects:
        base = Chem.RemoveHs(Chem.MolFromMolFile(f'{L.DOCK}/{c}_corr.sdf', removeHs=False))
        for tag, seed in (('conf1', 101), ('conf2', 202)):
            out = f'{L.REV}/ligands/conformers/{c}_{tag}'
            if not os.path.exists(out + '.pdbqt'):
                m = Chem.AddHs(Chem.MolFromSmiles(BL.COMPOUNDS[c][0])); p = AllChem.ETKDGv3(); p.randomSeed = seed
                assert AllChem.EmbedMolecule(m, p) == 0; AllChem.MMFFOptimizeMolecule(m)
                with Chem.SDWriter(out + '.sdf') as w: w.write(m)
                r = subprocess.run([mk, '-i', out + '.sdf', '-o', out + '.pdbqt'], capture_output=True, text=True); assert r.returncode == 0, r.stderr[:200]
            mm = Chem.RemoveHs(Chem.MolFromMolFile(out + '.sdf', removeHs=False))
            log.p(f'  conformer {c} {tag} (ETKDGv3 seed {seed}): best-aligned heavy-atom RMSD vs canonical {rdMolAlign.GetBestRMS(Chem.Mol(mm), Chem.Mol(base)):.2f} A')
    jobs = []
    for c in subjects:
        for t in ('fgfr3', 'fgfr1'):
            for s in SEEDS3:
                can = f'work/ligands/{c}_rebuilt.pdbqt'
                jobs.append(L.make_job('2d', can, V2[t], bc[t], box=30, exh=32, seed=s, cpu=0, conf='conf0', tag='2d exh32'))
                jobs.append(L.make_job('2d', can, V2[t], bc[t], box=20, exh=8, seed=s, cpu=0, conf='conf0', tag='2d box20'))
                jobs.append(L.make_job('2d', can, V2[t], bc[t], box=40, exh=8, seed=s, cpu=0, conf='conf0', tag='2d box40'))
                for tag in ('conf1', 'conf2'):
                    jobs.append(L.make_job('2d', f'work/ligands/conformers/{c}_{tag}.pdbqt', V2[t], bc[t], box=30, exh=8, seed=s, cpu=0, conf=tag, tag=f'2d {tag}'))
                if c in ALT: jobs.append(L.make_job('2d', f'work/ligands/alt_conformers/{c}_alt.pdbqt', V2[t], bc[t], box=30, exh=8, seed=s, cpu=0, conf='alt', tag='2d alt'))
    L.write_manifest(f'{L.REV}/manifests/d2.csv', jobs)
    json.dump(dict(mandatory=MANDATORY, boundary=boundary, subjects=subjects), open(f'{L.REV}/manifests/d2_subjects.json', 'w'))
    log.p(f'd2.csv: {len(jobs)} jobs for {len(subjects)} subjects (cost units ~ {sum(4 if j["exh"]==32 else 1 for j in jobs)})')
else:
    subj = json.load(open(f'{L.REV}/manifests/d2_subjects.json'))['subjects']
    meta = {}
    for m in L.all_meta():
        meta[(os.path.basename(m['ligand'])[:-6], os.path.basename(m['receptor'])[:-6], int(m['seed']), int(m['box']), int(m['exh']), m['id'].split('__')[-1])] = float(m['best_score'])
    summ = {r['compound']: r for r in csv.DictReader(open(f'{L.REV}/replicate_summary.csv'))}
    CONDS = [('base', 30, 8, 'conf0', '{c}_rebuilt'), ('exh32', 30, 32, 'conf0', '{c}_rebuilt'), ('box20', 20, 8, 'conf0', '{c}_rebuilt'), ('box40', 40, 8, 'conf0', '{c}_rebuilt'),
             ('conf1', 30, 8, 'conf1', '{c}_conf1'), ('conf2', 30, 8, 'conf2', '{c}_conf2'), ('alt', 30, 8, 'alt', '{c}_alt')]
    rows = []
    for c in subj:
        base_si = None; condsi = {}
        for name, box, exh, conf, lig in CONDS:
            if name == 'alt' and c not in ALT: continue
            vals = []
            for s in SEEDS3:
                a = meta.get((lig.format(c=c), 'receptor_fgfr3_v2', s, box, exh, conf)); b = meta.get((lig.format(c=c), 'receptor_fgfr1_v2', s, box, exh, conf))
                if a is None or b is None: vals = None; break
                vals.append((a, b, a - b))
            if vals is None: log.p(f'  incomplete {c} {name}'); continue
            m3 = sum(v[0] for v in vals) / 3; m1 = sum(v[1] for v in vals) / 3; si = [v[2] for v in vals]; msi = sum(si) / 3
            _, sdc = L.mean_sd(si)
            condsi[name] = msi
            if name == 'base': base_si = msi
            rows.append(dict(compound=c, condition=name, mean_fgfr3=round(m3, 3), mean_fgfr1=round(m1, 3), si_per_seed=' '.join(f'{x:+.2f}' for x in si), mean_si=round(msi, 3), sd_si_3seeds=round(sdc, 3)))
        sd5 = float(summ[c]['sd_si']); lo5, hi5 = float(summ[c]['ci_lo']), float(summ[c]['ci_hi']); m5 = float(summ[c]['mean_si'])
        rng = (max(condsi.values()) - min(condsi.values())) if condsi else float('nan')
        for r in rows:
            if r['compound'] == c:
                r['diff_from_base'] = round(r['mean_si'] - base_si, 3) if base_si is not None else ''
                r['range_over_conditions'] = round(rng, 3); r['sd_si_2c_five_seeds'] = sd5; r['outcome_2c'] = summ[c]['rule_outcome']
                r['sign_differs_from_2c_mean'] = 'yes' if (r['mean_si'] > 0) != (m5 > 0) else 'no'
                r['outside_2c_95_interval'] = 'yes' if (r['mean_si'] < lo5 or r['mean_si'] > hi5) else 'no'
                s3 = [float(x) for x in r['si_per_seed'].split()]; mm, sdd = L.mean_sd(s3); hw = 4.303 * sdd / math.sqrt(3)
                r['rule_adapted_3seed_t4.303'] = 'FGFR3-preferring' if mm + hw < 0 else ('FGFR1-preferring' if mm - hw > 0 else 'indistinguishable')
    with open(f'{L.REV}/sensitivity.csv', 'w', newline='') as f:
        keys = ['compound', 'condition', 'mean_fgfr3', 'mean_fgfr1', 'si_per_seed', 'mean_si', 'diff_from_base', 'range_over_conditions', 'sd_si_2c_five_seeds', 'outcome_2c', 'sign_differs_from_2c_mean', 'outside_2c_95_interval', 'rule_adapted_3seed_t4.303', 'sd_si_3seeds']
        w = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore'); w.writeheader(); w.writerows(rows)
    log.p(f'{"compound":16}{"cond":7}{"meanSI":>8}{"diff":>8}   (2c: mean {"":>6} SD)')
    for r in rows: log.p(f"{r['compound']:16}{r['condition']:7}{r['mean_si']:+8.3f}{r.get('diff_from_base', ''):>8}   sign-change {r.get('sign_differs_from_2c_mean')} outside-2c-CI {r.get('outside_2c_95_interval')} adapted-rule {r.get('rule_adapted_3seed_t4.303')}")
    for c in subj:
        rr = [r for r in rows if r['compound'] == c]
        if rr: log.p(f"  {c}: SI range over conditions {rr[0]['range_over_conditions']} (2c SD {rr[0]['sd_si_2c_five_seeds']}); conditions with sign change: {[r['condition'] for r in rr if r['sign_differs_from_2c_mean']=='yes']}; outside 2c interval: {[r['condition'] for r in rr if r['outside_2c_95_interval']=='yes']}")
