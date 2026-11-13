"""Rebuilds ligand SDF + PDBQT for all 21 compounds.
Canonical set: the SMILES strings in build_ligands.py -> ligands/<name>_rebuilt.*
Alternative conformers: embeddings from the isomeric SMILES in compounds.csv (different atom order) -> ligands/alt_conformers/<name>_alt.* (written for the compounds whose alternative differs).
Fixed conformer seed. Each rebuilt compound is compared with its existing *_corr.sdf / *_corr.pdbqt.
Pipeline (same settings as build_ligands.py): MolFromSmiles, AddHs, ETKDGv3(randomSeed=seed), MMFFOptimizeMolecule, SDWriter, mk_prepare_ligand.py (default options).
Usage: python phase1j_rebuild_ligands.py <docking_dir> <work_dir> [--seed 42]
Log: <work_dir>/ligands/rebuild_log.txt  Table: <work_dir>/ligands/rebuild_vs_corr.csv"""
import sys, os, csv, subprocess, platform, datetime, hashlib, importlib.util, shutil
import numpy as np
import rdkit, meeko
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolAlign

args = [a for a in sys.argv[1:] if not a.startswith('--')]
D, R = args[:2]
SEED = int(sys.argv[sys.argv.index('--seed') + 1]) if '--seed' in sys.argv else 42
if '--seed' in sys.argv: args = [a for a in args if a != str(SEED)]
os.chdir(D)
LD = f'{R}/ligands'; os.makedirs(LD, exist_ok=True)
log = open(f'{LD}/rebuild_log.txt', 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); log.write(s + '\n')
mk = shutil.which('mk_prepare_ligand.py')
P('# phase1j_rebuild_ligands.py', datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
P('# args:', sys.argv[1:], '| conformer seed', SEED)
P('# python', platform.python_version(), '| rdkit', rdBase.rdkitVersion, '| meeko', meeko.__version__, '| mk_prepare_ligand.py', mk)
P('# ETKDGv3 defaults except randomSeed; MMFF94 via MMFFOptimizeMolecule defaults (maxIters 200); Meeko CLI defaults')

spec = importlib.util.spec_from_file_location('build_ligands', f'{D}/build_ligands.py')
BL = importlib.util.module_from_spec(spec); spec.loader.exec_module(BL)

def embed(smiles):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = AllChem.ETKDGv3(); p.randomSeed = SEED
    cid = AllChem.EmbedMolecule(m, p)
    mm = AllChem.MMFFOptimizeMolecule(m)
    return m, cid, mm
def heavy(m): return Chem.RemoveHs(m)
def md5(p): return hashlib.md5(open(p, 'rb').read()).hexdigest()
def atom_lines(p): return [l for l in open(p) if l.startswith(('ATOM', 'HETATM'))]

rows = []
for r in csv.DictReader(open(f'{R}/compounds.csv')):
    n = r['compound']; smi = BL.COMPOUNDS[n][0]; smi_alt = r['isomeric_smiles']      # canonical = build_ligands.py string
    m, cid, mmff = embed(smi)
    sdf = f'{LD}/{n}_rebuilt.sdf'; pdbqt = f'{LD}/{n}_rebuilt.pdbqt'
    with Chem.SDWriter(sdf) as w: w.write(m)
    cp = subprocess.run([mk, '-i', sdf, '-o', pdbqt], capture_output=True, text=True)
    ok = cp.returncode == 0 and os.path.exists(pdbqt)
    # compare with existing _corr.sdf
    ex = Chem.MolFromMolFile(f'{n}_corr.sdf', removeHs=False); Chem.AssignStereochemistryFrom3D(ex)
    rb = Chem.MolFromMolFile(sdf, removeHs=False); Chem.AssignStereochemistryFrom3D(rb)
    key_ex, key_rb = Chem.MolToInchiKey(Chem.RemoveHs(ex)), Chem.MolToInchiKey(Chem.RemoveHs(rb))
    hx, hr = heavy(ex), heavy(rb)
    try: raw = rdMolAlign.CalcRMS(Chem.Mol(hr), Chem.Mol(hx))
    except Exception as e: raw = float('nan')
    try:
        pr = Chem.Mol(hr); best = rdMolAlign.GetBestRMS(pr, Chem.Mol(hx))
    except Exception as e: best = float('nan')
    # alternative: embedding from the compounds.csv canonical SMILES (different atom order)
    m2, cid2, _ = embed(smi_alt); m2h = heavy(Chem.RemoveHs(m2))
    try: raw2 = rdMolAlign.CalcRMS(Chem.Mol(m2h), Chem.Mol(hx))
    except Exception: raw2 = float('nan')
    try: best2 = rdMolAlign.GetBestRMS(Chem.Mol(m2h), Chem.Mol(hx))
    except Exception: best2 = float('nan')
    alt_written = ''
    if not (raw2 < 0.01):
        os.makedirs(f'{LD}/alt_conformers', exist_ok=True)
        asdf = f'{LD}/alt_conformers/{n}_alt.sdf'; apq = f'{LD}/alt_conformers/{n}_alt.pdbqt'
        with Chem.SDWriter(asdf) as w: w.write(m2)
        subprocess.run([mk, '-i', asdf, '-o', apq], capture_output=True, text=True); alt_written = os.path.relpath(apq, R)
    pq_same = ok and (md5(pdbqt) == md5(f'{n}_corr.pdbqt'))
    pq_atoms_same = ok and (atom_lines(pdbqt) == atom_lines(f'{n}_corr.pdbqt'))
    nrot = sum(1 for l in open(pdbqt) if l.startswith('BRANCH')) if ok else None
    rows.append(dict(compound=n, in_paper_20=r['in_paper_20'], rebuild_ok=ok, embed_cid=cid, mmff_return=mmff, inchikey_rebuilt=key_rb, inchikey_corr_sdf=key_ex,
                     inchikey_same='yes' if key_rb == key_ex else 'NO', heavy_rmsd_raw_noalign=round(raw, 3), heavy_rmsd_best_aligned=round(best, 3),
                     alt_conformer_from_compounds_csv_SMILES_raw_rmsd_vs_canonical=round(raw2, 3), alt_conformer_best_aligned_rmsd_vs_canonical=round(best2, 3), alt_conformer_pdbqt=alt_written, pdbqt_md5_identical=pq_same, pdbqt_atom_lines_identical=pq_atoms_same,
                     torsion_branches=nrot, smiles_used=smi, smiles_alt=smi_alt, seed=SEED, rdkit=rdBase.rdkitVersion, meeko=meeko.__version__,
                     rebuilt_sdf=os.path.relpath(sdf, R), rebuilt_pdbqt=os.path.relpath(pdbqt, R)))
    P(f"{n:16} ok={ok} InChIKey same={key_rb == key_ex}  heavy RMSD raw {raw:.3f}  best-aligned {best:.3f}  | alt (compounds.csv SMILES) raw {raw2:.3f} aligned {best2:.3f} {'WRITTEN '+alt_written if alt_written else ''} | pdbqt identical {pq_same}")
with open(f'{LD}/rebuild_vs_corr.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
diff = [x['compound'] for x in rows if x['inchikey_same'] != 'yes']
P('\ncompounds whose rebuilt InChIKey differs from the _corr.sdf InChIKey:', diff or 'none')
P('canonical rebuild vs existing _corr: max raw heavy-atom RMSD', max(x['heavy_rmsd_raw_noalign'] for x in rows))
P('alternative conformers written (differ from canonical):', [(x['compound'], x['alt_conformer_best_aligned_rmsd_vs_canonical']) for x in rows if x['alt_conformer_pdbqt']])
P('compounds whose rebuilt pdbqt atom lines are identical to the existing _corr.pdbqt:', sum(x['pdbqt_atom_lines_identical'] for x in rows), 'of', len(rows))
