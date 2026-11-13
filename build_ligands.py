"""build_ligands.py - builds the 21 docked compounds from SMILES.

Every compound is built on the 3,5-dichloropyridin-4-yl scaffold of dabogratinib, except custom1, which is an independent scaffold.
Preparation is identical for all: RDKit ETKDGv3 (randomSeed=42), MMFF94 optimisation, Meeko (mk_prepare_ligand.py).
Outputs <name>_corr.sdf and <name>_corr.pdbqt in the current folder.

Run: python build_ligands.py"""
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

# name: (SMILES_on_correct_scaffold, modification_in_original_terms, provenance)
COMPOUNDS = {
    # ── reference ────────────────────────────────────────────────────────────
    'dabogratinib':    ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CC5(C4)CN(S(C)(=O)=O)C5)nc3)c2c1)c1c(Cl)cncc1Cl',
                        'Reference: TYRA-300 (3,5-dichloropyridin-4-yl)', '9CD7 co-crystal ligand A1AV2'),
    # ── single-fragment variants ─────────────────────────────────────────────
    'v1_2F_pyridine':  ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CC5(C4)CN(S(C)(=O)=O)C5)nc3)c2c1)c1c(F)cncc1Cl',
                        '[A] 2-Cl->F (P-loop, smaller halogen)', 'single-fragment variant'),
    'v2_iPr_sulfonyl': ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CC5(C4)CN(S(C(C)C)(=O)=O)C5)nc3)c2c1)c1c(Cl)cncc1Cl',
                        '[E] Me->iPr sulfonamide (back pocket / gatekeeper)', 'single-fragment variant'),
    'v3_6F_indazole':  ('C[C\x40\x40H](Oc1cc(F)c2[nH]nc(-c3ccc(N4CC5(C4)CN(S(C)(=O)=O)C5)nc3)c2c1)c1c(Cl)cncc1Cl',
                        '[C] 6-F indazole (hinge electronics)', 'single-fragment variant'),
    'v4_Nme_indazole': ('C[C\x40\x40H](Oc1ccc2n(C)nc(-c3ccc(N4CC5(C4)CN(S(C)(=O)=O)C5)nc3)c2c1)c1c(Cl)cncc1Cl',
                        '[C] N1-methyl indazole (block hinge NH)', 'single-fragment variant'),
    'v5_pyrimidine':   ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3cnc(N4CC5(C4)CN(S(C)(=O)=O)C5)cn3)c2c1)c1c(Cl)cncc1Cl',
                        '[D] pyridine->pyrimidine linker (+HBA, polarity)', 'single-fragment variant'),
    # ── combinations built on the indazole-fluoro variant v3 ─────────────────
    'c3_v3v2':   ('C[C\x40\x40H](Oc1cc(F)c2[nH]nc(-c3ccc(N4CC5(C4)CN(S(C(C)C)(=O)=O)C5)nc3)c2c1)c1c(Cl)cncc1Cl',
                  '6-F indazole + iPr-sulfonyl', 'combination of variants'),
    'c3_v3v1':   ('C[C\x40\x40H](Oc1cc(F)c2[nH]nc(-c3ccc(N4CC5(C4)CN(S(C)(=O)=O)C5)nc3)c2c1)c1c(F)cncc1Cl',
                  '6-F indazole + 2-Cl->F', 'combination of variants'),
    'c3_v3v5':   ('C[C\x40\x40H](Oc1cc(F)c2[nH]nc(-c3cnc(N4CC5(C4)CN(S(C)(=O)=O)C5)cn3)c2c1)c1c(Cl)cncc1Cl',
                  '6-F indazole + pyrimidine linker', 'combination of variants'),
    'c3_v3v2v1': ('C[C\x40\x40H](Oc1cc(F)c2[nH]nc(-c3ccc(N4CC5(C4)CN(S(C(C)C)(=O)=O)C5)nc3)c2c1)c1c(F)cncc1Cl',
                  '6-F indazole + iPr-sulfonyl + 2-Cl->F', 'combination of variants'),
    'c3_v3v1v5': ('C[C\x40\x40H](Oc1cc(F)c2[nH]nc(-c3cnc(N4CC5(C4)CN(S(C)(=O)=O)C5)cn3)c2c1)c1c(F)cncc1Cl',
                  '6-F indazole + 2-Cl->F + pyrimidine', 'combination of variants'),
    # ── early candidates and hybrids ─────────────────────────────────────────
    'candidate1': ('CCS(=O)(=O)N1CC2(CN(c3ccc(-c4n[nH]c5ccc(O[C\x40H](C)c6c(Cl)cncc6Cl)cc45)cn3)C2)C1',
                   '[E] methylsulfonamide -> ethylsulfonamide', 'candidate1.py'),
    'candidate2': ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CC5(C4)CN(S(=O)(=O)C(F)(F)F)C5)nc3)c2c1)c1c(Cl)cncc1Cl',
                   '[E] methylsulfonyl -> trifluoromethylsulfonyl (CF3)', 'selectivity_pipeline.py'),
    'candidate3': ('CS(=O)(=O)N1CC2(CN(c3ccc(-c4n[nH]c5ccc(O[C\x40H](c6c(Cl)cncc6Cl)C6CC6)cc45)cn3)C2)C1',
                   '[B] methyl stereocentre -> cyclopropyl', 'selectivity_pipeline.py'),
    'custom1':    ('Cc1c(-c2cc(O[C\x40H](CO)c3ccc(F)cn3)c3c(Cl)cnn3c2)nnn1C1CCN(C#N)CC1',
                   'DE NOVO custom scaffold (NOT a TYRA-300 dichloropyridine analog): '
                   'CH2OH stereocentre, 2-fluoropyridine, chloro-pyrazolo core, '
                   'N1-methyltriazole, N-cyanamide piperidine', 'recovered from custom1.sdf (no script)'),
    'hybrid1':    ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CCC(C#N)CC4)nc3)c2c1)c1c(Cl)cncc1Cl',
                   '[E] diazaspiro[3.3]heptane-methylsulfonamide -> 4-cyanopiperidine',
                   'recovered from hybrid1.sdf (no script)'),
    'hybrid1_A':  ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CCC(C#N)(CO)CC4)nc3)c2c1)c1c(Cl)cncc1Cl',
                   '[E] 4-cyanopiperidine + 4-hydroxymethyl', 'recovered from hybrid1_A.sdf (no script)'),
    'hybrid1_B':  ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CC[C\x40H](C#N)[C\x40H](F)C4)nc3)c2c1)c1c(Cl)cncc1Cl',
                   '[E] 4-cyanopiperidine + 3-fluoro', 'recovered from hybrid1_B.sdf (no script)'),
    'hybrid1_C':  ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CCC(C#N)CC4)nc3)c2c1)c1c(F)cncc1Cl',
                   '[E] 4-cyanopiperidine + [A] 2-Cl->F pyridine', 'recovered from hybrid1_C.sdf (no script)'),
    'hybrid2_2F': ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CCC(C#N)CC4)nc3)c2c1)c1c(F)cncc1F',
                   '[E] 4-cyanopiperidine + [A] 2,3-diCl->diF pyridine', 'recovered from hybrid2_2F.sdf (no script)'),
    'hybrid2_OH': ('C[C\x40\x40H](Oc1ccc2[nH]nc(-c3ccc(N4CC[C\x40\x40H](O)[C\x40H](C#N)CC4)nc3)c2c1)c1c(F)cncc1Cl',
                   '[E] 4-cyanopiperidine + 3-hydroxy + [A] 2-Cl->F pyridine', 'recovered from hybrid2_OH.sdf (no script)'),
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
