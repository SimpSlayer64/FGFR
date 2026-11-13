"""Receptor re-preparation test: prepares receptors into work/prep/.
noH     : obabel <chainA.pdb> -O <out> -xrn
rep1..3 : obabel <chainA.pdb> -O <out> -xr -p 7.4   (as in build_receptors.py), each re-run from the CIF
Writes manifests/g2.csv. Resumable (existing outputs are kept). Usage: PYTHONDONTWRITEBYTECODE=1 python p2g_prep.py"""
import os, sys, json, subprocess, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_lib as L
import gemmi, numpy as np
log = L.Log('p2g_prep')
bc = L.box_centres(); LIG = 'A1AV2'
def extract(cif, chain, out_pdb):
    st = gemmi.read_structure(f'{L.DOCK}/{cif}'); st.setup_entities(); model = st[0]
    prot = gemmi.Structure(); prot.add_model(gemmi.Model('1')); ch = gemmi.Chain(chain); lig = []
    for c in model:
        if c.name != chain: continue
        for res in c:
            info = gemmi.find_tabulated_residue(res.name)
            if res.name == LIG: lig += [[a.pos.x, a.pos.y, a.pos.z] for a in res]
            elif info and info.is_amino_acid(): ch.add_residue(res)
    prot[0].add_chain(ch); prot.write_pdb(out_pdb)
    c = np.array(lig).mean(0); return [round(float(c[0]), 2), round(float(c[1]), 2), round(float(c[2]), 2)]
ARMS = {'noH': ['-xrn'], 'rep1': ['-xr', '-p', '7.4'], 'rep2': ['-xr', '-p', '7.4'], 'rep3': ['-xr', '-p', '7.4']}
jobs = []
for arm, flags in ARMS.items():
    d = f'{L.REV}/prep/{arm}'; os.makedirs(d, exist_ok=True)
    for tag, cif in (('fgfr3', '9CD7.cif'), ('fgfr1', '9CD5.cif')):
        pdb = f'{d}/{tag}_chainA.pdb'; out = f'{d}/receptor_{tag}_{arm}.pdbqt'
        centre = extract(cif, 'A', pdb)
        assert centre == bc[tag], (centre, bc[tag])
        same = hashlib.md5(open(pdb, 'rb').read()).hexdigest() == hashlib.md5(open(f'{L.DOCK}/{tag}_chainA.pdb', 'rb').read()).hexdigest()
        if not os.path.exists(out):
            cmd = ['obabel', pdb, '-O', out] + flags
            r = subprocess.run(cmd, capture_output=True, text=True); log.p(' '.join(cmd), '->', (r.stdout + r.stderr).strip().split('\n')[-1])
            open(f'{d}/{tag}_obabel.log', 'w').write(' '.join(cmd) + '\n' + r.stdout + r.stderr)
        n = sum(1 for l in open(out) if l.startswith(('ATOM', 'HETATM'))); nh = sum(1 for l in open(out) if l.startswith(('ATOM', 'HETATM')) and l[77:79].strip() == 'HD')
        log.p(f'{arm} {tag}: chain-A pdb identical to project {tag}_chainA.pdb: {same}; centre {centre} == box_centers.json; pdbqt {n} atoms, {nh} HD')
        for s in L.SEEDS5:
            jobs.append(L.make_job('prep', 'work/ligands/dabogratinib_rebuilt.pdbqt', os.path.relpath(out, L.DOCK), bc[tag], seed=s, cpu=0, conf='conf0', tag=f'2g {arm} {tag}'))
L.write_manifest(f'{L.REV}/manifests/g2.csv', jobs); log.p(f'g2.csv: {len(jobs)} jobs')
# how different are the three re-preps (polar-H placement)?
def hd(p): return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in open(p) if l.startswith('ATOM') and l[77:79].strip() == 'HD'])
for tag in ('fgfr3', 'fgfr1'):
    A = [hd(f'{L.REV}/prep/rep{i}/receptor_{tag}_rep{i}.pdbqt') for i in (1, 2, 3)] + [hd(f'{L.DOCK}/receptor_{tag}_v2.pdbqt')]
    names = ['rep1', 'rep2', 'rep3', 'v2']
    for i in range(4):
        for j in range(i + 1, 4):
            if A[i].shape == A[j].shape:
                d = np.linalg.norm(A[i] - A[j], axis=1); log.p(f'  {tag} polar-H positions {names[i]} vs {names[j]}: {len(d)} HD atoms, moved > 0.01 A: {int((d>0.01).sum())}, max shift {d.max():.2f} A, mean {d.mean():.3f} A')
            else: log.p(f'  {tag} {names[i]} vs {names[j]}: different HD counts {A[i].shape[0]} vs {A[j].shape[0]}')
