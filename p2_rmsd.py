"""Symmetry-aware heavy-atom RMSD of docked poses against every crystal copy of a ligand.
One script; the CIF, chain and ligand residue are arguments.

  python p2_rmsd.py --cif 9CD7.cif --chain A --ligand-resn A1AV2 --template-smiles "<SMILES>" --poses a.pdbqt b.pdbqt ... --out rmsd_part.csv [--core-ref]

Method: crystal copies = every chain's residue --ligand-resn in the CIF; bond orders assigned from --template-smiles (AssignBondOrdersFromTemplate);
copies in chains other than --chain are moved into the --chain frame by superposing that chain's protein CA onto --chain's (identical residue number and name,
iterative 2.0 A trimming), so all RMSDs are computed in the docking receptor frame. Poses are read with Meeko (RDKitMolCreate.from_pdbqt_mol), hydrogens removed.
RMSD = rdkit.Chem.rdMolAlign.CalcRMS(probe, ref, ...): computed in place, no alignment, symmetry-aware. If the pose is the same molecule as the crystal ligand the
whole heavy-atom molecule is compared ('full'). With --core-ref (analogs) the atoms compared are the maximum common substructure between pose and crystal ligand
(rdFMCS: elements, any bond, rings match rings; heavy atoms), identical rule for every receptor; all symmetry-equivalent mappings are passed as `map`.
A core of fewer than --min-core heavy atoms (default 20) is reported as 'core_too_small' and no RMSD is given.
"""
import os, sys, csv, argparse, itertools, datetime, platform
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, gemmi
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolAlign, rdFMCS, AllChem
from meeko import PDBQTMolecule, RDKitMolCreate
import meeko

ap = argparse.ArgumentParser()
ap.add_argument('--cif', required=True); ap.add_argument('--chain', required=True); ap.add_argument('--ligand-resn', required=True)
ap.add_argument('--template-smiles', required=True); ap.add_argument('--poses', nargs='+', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--core-ref', action='store_true'); ap.add_argument('--min-core', type=int, default=20); ap.add_argument('--max-maps', type=int, default=5000)
a = ap.parse_args()
print(f'# p2_rmsd.py {datetime.datetime.now():%Y-%m-%d %H:%M:%S} args={sys.argv[1:]}\n# python {platform.python_version()} rdkit {rdBase.rdkitVersion} meeko {meeko.__version__} gemmi {gemmi.__version__}', flush=True)

# ---------------- crystal copies in the docking-chain frame
st = gemmi.read_structure(a.cif); st.setup_entities()
def prot(ch): return {r.seqid.num: r for r in ch if gemmi.find_tabulated_residue(r.name) and gemmi.find_tabulated_residue(r.name).is_amino_acid()}
chains = {c.name: c for c in st[0]}
P0 = prot(chains[a.chain])
def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0); H = np.matmul((Q - qc).T, P - pc); U, S, Vt = np.linalg.svd(H); d = np.sign(np.linalg.det(np.matmul(Vt.T, U.T)))
    R = np.matmul(np.matmul(Vt.T, np.diag([1, 1, d])), U.T); return R, pc - np.matmul(R, qc)
def frame(ch):
    if ch == a.chain: return np.eye(3), np.zeros(3), len(P0), 0.0
    Px = prot(chains[ch]); pr = [n for n in P0 if n in Px and P0[n].name == Px[n].name]
    for _ in range(50):
        A = np.array([P0[n]['CA'][0].pos.tolist() for n in pr]); B = np.array([Px[n]['CA'][0].pos.tolist() for n in pr])
        R, t = kabsch(A, B); d = np.linalg.norm(A - (np.matmul(B, R.T) + t), axis=1); keep = [n for n, di in zip(pr, d) if di <= 2.0]
        if len(keep) == len(pr): break
        pr = keep
    return R, t, len(pr), float(np.sqrt((d ** 2).mean()))
template = Chem.MolFromSmiles(a.template_smiles)
copies = {}
for cname, ch in chains.items():
    for res in ch:
        if res.name != a.ligand_resn: continue
        R, t, n, rms = frame(cname)
        xyz = np.matmul(np.array([x.pos.tolist() for x in res]), R.T) + t
        lines = [f'HETATM{i+1:5d} {x.name:<4s} LIG A   1    {xyz[i][0]:8.3f}{xyz[i][1]:8.3f}{xyz[i][2]:8.3f}  1.00  0.00          {x.element.name:>2s}' for i, x in enumerate(res)]
        m = Chem.MolFromPDBBlock('\n'.join(lines), sanitize=False, removeHs=False, proximityBonding=True)
        m = AllChem.AssignBondOrdersFromTemplate(template, m); Chem.SanitizeMol(m)
        copies[cname] = Chem.RemoveHs(m)
        print(f'crystal copy chain {cname}: {m.GetNumAtoms()} heavy atoms, chain->{a.chain} frame from {n} CA pairs (RMSD {rms:.2f} A)', flush=True)
assert copies, 'no ligand copies found'

def pose_mol(path):
    pm = PDBQTMolecule.from_file(path, skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(pm)
    m = mols[0]; n = m.GetNumConformers()
    return Chem.RemoveHs(m), n
def centroid(m, cid): return m.GetConformer(cid).GetPositions().mean(0)

rows = []
core_cache = {}
for pf in a.poses:
    try: pm, nconf = pose_mol(pf)
    except Exception as e:
        print('POSE READ FAILED', pf, type(e).__name__, e, flush=True); rows.append(dict(pose_file=pf, model=0, status=f'read_failed {e}')); continue
    for cname, ref in copies.items():
        key = (pm.GetNumAtoms(), Chem.MolToSmiles(pm), cname)
        if a.core_ref:
            if key not in core_cache:
                res = rdFMCS.FindMCS([pm, ref], atomCompare=rdFMCS.AtomCompare.CompareElements, bondCompare=rdFMCS.BondCompare.CompareAny,
                                     ringMatchesRingOnly=True, completeRingsOnly=False, timeout=20)
                q = Chem.MolFromSmarts(res.smartsString)
                mp = pm.GetSubstructMatches(q, uniquify=False, maxMatches=2000); mr = ref.GetSubstructMatches(q, uniquify=False, maxMatches=2000)
                maps = [[(p, r) for p, r in zip(x, y)] for x in mp for y in mr][:a.max_maps] if res.numAtoms >= a.min_core and mp and mr else []
                core_cache[key] = (res.numAtoms, maps, len(mp), len(mr))
            ncore, maps, nmp, nmr = core_cache[key]
        for cid in range(nconf):
            row = dict(pose_file=pf, model=cid + 1, crystal_copy=cname, status='ok', n_pose_heavy=pm.GetNumAtoms())
            d = float(np.linalg.norm(centroid(pm, cid) - centroid(ref, 0))); row['centroid_dist_heavy'] = round(d, 3)
            try:
                if a.core_ref:
                    row['core_atoms'] = ncore; row['n_maps'] = len(maps)
                    if not maps: row['status'] = 'core_too_small'; row['rmsd'] = ''
                    else: row['rmsd'] = round(rdMolAlign.CalcRMS(pm, ref, cid, 0, maps), 3)
                else:
                    row['core_atoms'] = pm.GetNumAtoms(); row['n_maps'] = ''
                    row['rmsd'] = round(rdMolAlign.CalcRMS(pm, ref, cid, 0), 3)
            except Exception as e:
                row['status'] = f'calcrms_failed {type(e).__name__}: {e}'; row['rmsd'] = ''
            rows.append(row)
    print(f'{os.path.basename(pf)}: {nconf} poses done', flush=True)
keys = ['pose_file', 'model', 'crystal_copy', 'status', 'rmsd', 'core_atoms', 'n_maps', 'n_pose_heavy', 'centroid_dist_heavy']
with open(a.out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore'); w.writeheader(); w.writerows(rows)
bad = [r for r in rows if str(r.get('status', '')).startswith('calcrms_failed') or str(r.get('status', '')).startswith('read_failed')]
print(f'wrote {a.out}: {len(rows)} rows, failures {len(bad)}', flush=True)
if bad: sys.exit(3)
