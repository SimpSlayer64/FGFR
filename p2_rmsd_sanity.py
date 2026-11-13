"""Sanity checks of the RMSD method: self-RMSD is 0; a pose with the two symmetry-equivalent chlorine-bearing pyridine carbons
(and their Cl and CH neighbours) swapped gives about 0 with CalcRMS (symmetry-aware) but a non-zero value with a plain atom-order RMSD."""
import os, sys, numpy as np, gemmi
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolAlign
D = os.environ.get('FGFR3_DOCK', os.getcwd())
print('# p2_rmsd_sanity', 'rdkit', rdBase.rdkitVersion)
tmpl = Chem.MolFromSmiles('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CC5(C4)CN(S(C)(=O)=O)C5)nc3)c2c1)c1c(Cl)cncc1Cl')
st = gemmi.read_structure(f'{D}/9CD7.cif'); res = [r for r in st[0]['A'] if r.name == 'A1AV2'][0]
xyz = np.array([x.pos.tolist() for x in res])
lines = [f'HETATM{i+1:5d} {x.name:<4s} LIG A   1    {xyz[i][0]:8.3f}{xyz[i][1]:8.3f}{xyz[i][2]:8.3f}  1.00  0.00          {x.element.name:>2s}' for i, x in enumerate(res)]
m = Chem.MolFromPDBBlock('\n'.join(lines), sanitize=False, proximityBonding=True); m = AllChem.AssignBondOrdersFromTemplate(tmpl, m); Chem.SanitizeMol(m); ref = Chem.RemoveHs(m)
print('self RMSD (CalcRMS):', round(rdMolAlign.CalcRMS(Chem.Mol(ref), ref), 6))
# flip: reflect the dichloropyridyl ring about the C4-C(ring axis): swap coordinates of the two symmetry-equivalent ortho carbons, their Cl, and the two meta CH
patt = Chem.MolFromSmarts('[c;$(c(Cl)c(c)C)](Cl)')
ar = [a for a in ref.GetAtoms() if a.GetSymbol() == 'Cl']
c1, c2 = [n.GetIdx() for cl in ar for n in cl.GetNeighbors()]
cl1, cl2 = [a.GetIdx() for a in ar]
ring_c = [n.GetIdx() for n in ref.GetAtomWithIdx(c1).GetNeighbors() if n.GetIdx() != cl1 and n.GetIsAromatic()]
ring_c2 = [n.GetIdx() for n in ref.GetAtomWithIdx(c2).GetNeighbors() if n.GetIdx() != cl2 and n.GetIsAromatic()]
# meta CH neighbours (not the attachment carbon): neighbours of c1/c2 that are bonded to the ring N side
pos = ref.GetConformer().GetPositions().copy()
def swap(i, j): pos[[i, j]] = pos[[j, i]]
swap(c1, c2); swap(cl1, cl2)
meta = [n for n in ring_c if ref.GetAtomWithIdx(n).GetDegree() == 2]; meta2 = [n for n in ring_c2 if ref.GetAtomWithIdx(n).GetDegree() == 2]
if meta and meta2: swap(meta[0], meta2[0])
flip = Chem.Mol(ref); cf = flip.GetConformer()
for i in range(flip.GetNumAtoms()): cf.SetAtomPosition(i, pos[i].tolist())
plain = float(np.sqrt(((pos - ref.GetConformer().GetPositions()) ** 2).sum(1).mean()))
print('flipped-pyridine pose: plain atom-order RMSD', round(plain, 3), '| CalcRMS (symmetry-aware)', round(rdMolAlign.CalcRMS(flip, ref), 6))
