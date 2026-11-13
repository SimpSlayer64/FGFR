"""build_ligands.py - builds the 21 docked compounds from SMILES.

The compounds are dabogratinib and 20 analogs defined by SMILES; custom1 is an independent scaffold and the others are modifications of dabogratinib.
Preparation is identical for all: RDKit ETKDGv3 (randomSeed=42), MMFF94 optimisation, Meeko (mk_prepare_ligand.py).
Outputs <name>_corr.sdf and <name>_corr.pdbqt in the current folder.

Run: python build_ligands.py"""
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

# name: (SMILES, chemical change relative to dabogratinib, provenance)
COMPOUNDS = {
    # parent
    'dabogratinib': ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CC5(C4)CN(S(C)(=O)=O)C5)nc3)c2c1)c1c(Cl)cncc1Cl',
        'Parent compound: 3,5-dichloropyridin-4-yl, methylsulfonyl-diazaspiroheptane',
        '9CD7 co-crystal ligand A1AV2'),
    # single changes to dabogratinib
    'v1_2F_pyridine': ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CC5(C4)CN(S(C)(=O)=O)C5)nc3)c2c1)c1c(F)cncc1Cl',
        'one Cl to F on the pyridine (3-chloro-5-fluoropyridin-4-yl)',
        'defined in an earlier script (not shipped)'),
    'v2_iPr_sulfonyl': ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CC5(C4)CN(S(C(C)C)(=O)=O)C5)nc3)c2c1)c1c(Cl)cncc1Cl',
        'methylsulfonyl to isopropylsulfonyl',
        'defined in an earlier script (not shipped)'),
    'v3_6F_indazole': ('C[C\x40\x40H](Oc1cc(F)c2[nH]nc(-c3ccc(N4CC5(C4)CN(S(C)(=O)=O)C5)nc3)c2c1)c1c(Cl)cncc1Cl',
        'F added on the indazole benzo ring, at position 7',
        'defined in an earlier script (not shipped)'),
    'v4_Nme_indazole': ('C[C\x40\x40H](Oc1ccc2n(C)nc(-c3ccc(N4CC5(C4)CN(S(C)(=O)=O)C5)nc3)c2c1)c1c(Cl)cncc1Cl',
        'methyl on indazole N1',
        'defined in an earlier script (not shipped)'),
    'v5_pyrimidine': ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3cnc(N4CC5(C4)CN(S(C)(=O)=O)C5)cn3)c2c1)c1c(Cl)cncc1Cl',
        'pyridine ring between indazole and spiro amine replaced by pyrazine',
        'defined in an earlier script (not shipped)'),
    # combinations built on the 7-F indazole variant
    'c3_v3v2': ('C[C\x40\x40H](Oc1cc(F)c2[nH]nc(-c3ccc(N4CC5(C4)CN(S(C(C)C)(=O)=O)C5)nc3)c2c1)c1c(Cl)cncc1Cl',
        '7-F indazole + methylsulfonyl to isopropylsulfonyl',
        'defined in an earlier script (not shipped)'),
    'c3_v3v1': ('C[C\x40\x40H](Oc1cc(F)c2[nH]nc(-c3ccc(N4CC5(C4)CN(S(C)(=O)=O)C5)nc3)c2c1)c1c(F)cncc1Cl',
        '7-F indazole + one Cl to F on the pyridine',
        'defined in an earlier script (not shipped)'),
    'c3_v3v5': ('C[C\x40\x40H](Oc1cc(F)c2[nH]nc(-c3cnc(N4CC5(C4)CN(S(C)(=O)=O)C5)cn3)c2c1)c1c(Cl)cncc1Cl',
        '7-F indazole + pyrazine in place of the pyridine ring',
        'defined in an earlier script (not shipped)'),
    'c3_v3v2v1': ('C[C\x40\x40H](Oc1cc(F)c2[nH]nc(-c3ccc(N4CC5(C4)CN(S(C(C)C)(=O)=O)C5)nc3)c2c1)c1c(F)cncc1Cl',
        '7-F indazole + isopropylsulfonyl + one Cl to F on the pyridine',
        'defined in an earlier script (not shipped)'),
    'c3_v3v1v5': ('C[C\x40\x40H](Oc1cc(F)c2[nH]nc(-c3cnc(N4CC5(C4)CN(S(C)(=O)=O)C5)cn3)c2c1)c1c(F)cncc1Cl',
        '7-F indazole + one Cl to F on the pyridine + pyrazine in place of the pyridine ring',
        'defined in an earlier script (not shipped)'),
    # other analogs and hybrids
    'candidate1': ('CCS(=O)(=O)N1CC2(CN(c3ccc(-c4n[nH]c5ccc(O[C\x40H](C)c6c(Cl)cncc6Cl)cc45)cn3)C2)C1',
        'methylsulfonyl to ethylsulfonyl',
        'defined in candidate1.py (not shipped)'),
    'candidate2': ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CC5(C4)CN(S(=O)(=O)C(F)(F)F)C5)nc3)c2c1)c1c(Cl)cncc1Cl',
        'methylsulfonyl to trifluoromethylsulfonyl',
        'defined in selectivity_pipeline.py (not shipped)'),
    'candidate3': ('CS(=O)(=O)N1CC2(CN(c3ccc(-c4n[nH]c5ccc(O[C\x40H](c6c(Cl)cncc6Cl)C6CC6)cc45)cn3)C2)C1',
        'methyl on the ether carbon to cyclopropyl',
        'defined in selectivity_pipeline.py (not shipped)'),
    'custom1': ('Cc1c(-c2cc(O[C\x40H](CO)c3ccc(F)cn3)c3c(Cl)cnn3c2)nnn1C1CCN(C#N)CC1',
        'Independent scaffold: chloro-pyrazolo[1,5-a]pyridine core, ether to a 1-(5-fluoropyridin-2-yl)-2-hydroxyethyl group, methyl-1,2,3-triazole, N-cyano piperidine',
        'structure read from the existing SDF file custom1.sdf; no generating script exists'),
    'hybrid1': ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CCC(C#N)CC4)nc3)c2c1)c1c(Cl)cncc1Cl',
        'methylsulfonyl-diazaspiroheptane to 4-cyanopiperidine',
        'structure read from the existing SDF file hybrid1.sdf; no generating script exists'),
    'hybrid1_A': ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CCC(C#N)(CO)CC4)nc3)c2c1)c1c(Cl)cncc1Cl',
        '4-cyanopiperidine with 4-hydroxymethyl',
        'structure read from the existing SDF file hybrid1_A.sdf; no generating script exists'),
    'hybrid1_B': ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CC[C\x40H](C#N)[C\x40H](F)C4)nc3)c2c1)c1c(Cl)cncc1Cl',
        '4-cyanopiperidine with 3-fluoro',
        'structure read from the existing SDF file hybrid1_B.sdf; no generating script exists'),
    'hybrid1_C': ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CCC(C#N)CC4)nc3)c2c1)c1c(F)cncc1Cl',
        '4-cyanopiperidine + one Cl to F on the pyridine',
        'structure read from the existing SDF file hybrid1_C.sdf; no generating script exists'),
    'hybrid2_2F': ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CCC(C#N)CC4)nc3)c2c1)c1c(F)cncc1F',
        '4-cyanopiperidine + both Cl to F on the pyridine (3,5-difluoropyridin-4-yl)',
        'structure read from the existing SDF file hybrid2_2F.sdf; no generating script exists'),
    'hybrid2_OH': ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CC[C\x40\x40H](O)[C\x40H](C#N)CC4)nc3)c2c1)c1c(F)cncc1Cl',
        'seven-membered azepane ring with hydroxy and cyano substituents + one Cl to F on the pyridine',
        'structure read from the existing SDF file hybrid2_OH.sdf; no generating script exists'),
}


def lipinski(m):
    mw, logp = Descriptors.MolWt(m), Descriptors.MolLogP(m)
    hbd, hba = Descriptors.NumHDonors(m), Descriptors.NumHAcceptors(m)
    f = []
    if mw > 500: f.append(f'MW={mw:.0f}')
    if logp > 5: f.append(f'logP={logp:.1f}')
    if hbd > 5:  f.append(f'HBD={hbd}')
    if hba > 10: f.append(f'HBA={hba}')
    return f


def build(name, smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f'invalid SMILES for {name}')
    mh = Chem.AddHs(m)
    params = AllChem.ETKDGv3(); params.randomSeed = 42
    if AllChem.EmbedMolecule(mh, params) == -1:
        raise RuntimeError(f'embedding failed for {name}')
    AllChem.MMFFOptimizeMolecule(mh)
    sdf = f'{name}_corr.sdf'
    with Chem.SDWriter(sdf) as w:
        w.write(mh)
    pdbqt = f'{name}_corr.pdbqt'
    r = subprocess.run(['mk_prepare_ligand.py', '-i', sdf, '-o', pdbqt],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f'meeko failed for {name}: {r.stderr[:200]}')
    return m


if __name__ == '__main__':
    for name, (smiles, mod, prov) in COMPOUNDS.items():
        m = build(name, smiles)
        flags = lipinski(m)
        print(f'{name:14} {rdMolDescriptors.CalcMolFormula(m):22} '
              f'{"; ".join(flags) if flags else "Lipinski OK":22} | {mod}  [{prov}]')
    # sanity: corrected dabogratinib must equal crystal TYRA-300
    dabo = Chem.MolFromSmiles(COMPOUNDS['dabogratinib'][0])
    assert Chem.MolToInchiKey(dabo) == 'JOAFWIHZEBKYQK-OAHLLOKOSA-N', 'dabogratinib != TYRA-300!'
    print(f'\nBuilt {len(COMPOUNDS)} ligands. dabogratinib == TYRA-300 (verified).')
