"""build_receptors.py - builds the FGFR3 and FGFR1 receptor pair.

  FGFR3 from 9CD7 chain A, FGFR1 from 9CD5 chain A (both with the same inhibitor bound).
Both are prepared by the same steps: keep the chain-A protein residues, drop ligand, water and ions, convert with
  obabel chain.pdb -O receptor.pdbqt -xr -p 7.4
and set the box centre to the centroid of the chain-A ligand (residue A1AV2) of the same structure.
Reads 9CD7.cif and 9CD5.cif from the current folder; writes receptor_fgfr3_v2.pdbqt, receptor_fgfr1_v2.pdbqt and box_centers.json.

Open Babel places polar hydrogens in a way that is not bitwise reproducible, so a rebuilt receptor is equivalent to, but not identical with,
the receptor files supplied in data/. Use the supplied files to reproduce the stored scores.

Run: python build_receptors.py"""
import subprocess, json, os
import gemmi
import numpy as np

LIG = 'A1AV2'
JOBS = [('fgfr3', '9CD7.cif', 'A'),
        ('fgfr1', '9CD5.cif', 'A')]

centers = {}
for tag, cif, chain in JOBS:
    st = gemmi.read_structure(cif)
    st.setup_entities()
    model = st[0]
    prot = gemmi.Structure(); prot.add_model(gemmi.Model('1'))
    ch = gemmi.Chain(chain)
    lig_xyz = []
    for c in model:
        if c.name != chain:
            continue
        for res in c:
            info = gemmi.find_tabulated_residue(res.name)
            if res.name == LIG:
                lig_xyz += [[a.pos.x, a.pos.y, a.pos.z] for a in res]
            elif info and info.is_amino_acid():
                ch.add_residue(res)
    prot[0].add_chain(ch)
    prot_pdb = f'{tag}_chainA.pdb'
    prot.write_pdb(prot_pdb)

    c = np.array(lig_xyz).mean(0)
    centers[tag] = [round(float(c[0]), 2), round(float(c[1]), 2), round(float(c[2]), 2)]

    out = f'receptor_{tag}_v2.pdbqt'
    subprocess.run(['obabel', prot_pdb, '-O', out, '-xr', '-p', '7.4'],
                   capture_output=True, text=True)
    natoms = sum(1 for l in open(out) if l.startswith(('ATOM', 'HETATM'))) if os.path.exists(out) else 0
    print(f'{tag}: chainA residues={len(ch)}  -> {out} ({natoms} atoms)  box center={centers[tag]}')

json.dump(centers, open('box_centers.json', 'w'), indent=1)
print('box centers written to box_centers.json:', centers)
