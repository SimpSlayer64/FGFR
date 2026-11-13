"""Builds the set of unrelated ATP-site ligands (no docking): the ATP-site ligand of every ensemble entry (phase3b/receptor_list.csv).
Excluded: dabogratinib, covalent ligands, metal-containing ligands, molecular weight > 600 Da.
Ligands are built from the component-dictionary SMILES (RCSB data API; retrieval time recorded in raw/): RDKit MolFromSmiles, AddHs, ETKDGv3 randomSeed 42, MMFF optimisation, SDF, Meeko mk_prepare_ligand.py defaults.
Covalent = (i) a covale record in the entry's _struct_conn joining the ligand to a protein residue in any copy, or (ii) a reactive-warhead SMARTS match on the CCD SMILES (acrylamide, vinyl sulfonamide, alpha-halo carbonyl, propiolamide, nitrile-activated acrylamide), listed for inspection.
Usage: PYTHONDONTWRITEBYTECODE=1 python p3c_null_build.py"""
import sys, os, csv, json, time, shutil, subprocess, hashlib, datetime, urllib.request, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_lib as L
import gemmi
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
log = L.Log('p3c_null_build'); C = f'{L.REV}/phase3c'; B = f'{L.REV}/phase3b'; P3 = f'{L.REV}/phase3'
MAXDA = 600.0
NONMETAL = {'H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'Te', 'I', 'D'}
WARHEADS = {'acrylamide/enone (C=C-C=O, terminal or substituted)': '[CX3;!R]=[CX3;!R]-[CX3](=O)[#7,#6]', 'vinyl sulfonamide/sulfone': '[CX3]=[CX3]-[SX4](=O)(=O)', 'alpha-halo carbonyl': '[CX4]([Cl,Br,I])[CX3]=O',
            'propiolamide/ynone': '[CX2]#[CX2]-[CX3]=O', 'cyanoacrylamide': 'N#C-C(=C)-C=O', 'acyl/aryl fluorosulfonate or sulfonyl fluoride': '[SX4](=O)(=O)F', 'epoxide/aziridine': '[C;r3]1[O,N][C;r3]1', 'isothiocyanate': 'N=C=S'}
def get(url, path):
    if not os.path.exists(path):
        for k in range(3):
            try:
                d = urllib.request.urlopen(url, timeout=60).read(); open(path, 'wb').write(d); json.dump(dict(url=url, retrieved=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')), open(path + '.meta', 'w')); break
            except Exception as ex: last = ex; time.sleep(2)
        else: raise last
    return json.load(open(path))
RL = list(csv.DictReader(open(f'{B}/receptor_list.csv')))
lig_entries = {}
for r in RL: lig_entries.setdefault(r['ligand_resname'], []).append(r)
log.p(f'{len(RL)} ensemble entries, {len(lig_entries)} distinct ATP-site ligand IDs')
def covalent_records(entry, comp):
    b = gemmi.cif.read(f'{P3}/raw/pdb/{entry}.cif').sole_block(); t = b.find('_struct_conn.', ['conn_type_id', 'ptnr1_label_comp_id', 'ptnr2_label_comp_id', 'ptnr1_auth_asym_id', 'ptnr2_auth_asym_id', 'ptnr1_auth_seq_id', 'ptnr2_auth_seq_id', 'ptnr1_label_atom_id', 'ptnr2_label_atom_id'])
    out = []
    for r in t:
        if r[0] != 'covale': continue
        a, c = r[1], r[2]
        if (a == comp) != (c == comp):
            other = c if a == comp else a
            if gemmi.find_tabulated_residue(other) and gemmi.find_tabulated_residue(other).is_amino_acid(): out.append(f'{r[3]}:{r[1]}{r[5]}.{r[7]}-{r[4]}:{r[2]}{r[6]}.{r[8]}')
    return out
dab = Chem.MolFromSmiles([r for r in csv.DictReader(open(f'{L.REV}/compounds.csv')) if r['compound'] == 'dabogratinib'][0]['isomeric_smiles'])
dab_ik = Chem.MolToInchiKey(dab).split('-')[0]; dab_fp = AllChem.GetMorganFingerprintAsBitVect(dab, 2, 2048)
mk = shutil.which('mk_prepare_ligand.py'); log.p('mk_prepare_ligand.py', mk, '| dabogratinib InChIKey block 1', dab_ik)
rows = []
for comp, ents in sorted(lig_entries.items()):
    d = get(f'https://data.rcsb.org/rest/v1/core/chemcomp/{comp}', f'{C}/raw/ccd_{comp}.json'); cc = d['chem_comp']; desc = d.get('rcsb_chem_comp_descriptor', {})
    smi = desc.get('SMILES_stereo') or desc.get('SMILES'); mw = float(cc['formula_weight']); formula = cc['formula']; elems = set(re.findall(r'[A-Z][a-z]?', formula))
    metals = sorted(elems - NONMETAL); cov = {e['entry']: covalent_records(e['entry'], comp) for e in ents}; cov = {k: v for k, v in cov.items() if v}
    m = Chem.MolFromSmiles(smi); heavy = m.GetNumHeavyAtoms() if m else ''
    wh = [n for n, s in WARHEADS.items() if m and m.HasSubstructMatch(Chem.MolFromSmarts(s))]
    ik = Chem.MolToInchiKey(m).split('-')[0] if m else ''
    reasons = []
    if comp == 'A1AV2' or ik == dab_ik: reasons.append('TYRA-300 / dabogratinib')
    if cov: reasons.append('covalent (struct_conn): ' + str(cov))
    if wh: reasons.append('covalent warhead SMARTS: ' + '; '.join(wh))
    if metals: reasons.append('metal-containing: ' + ','.join(metals))
    if mw > MAXDA: reasons.append(f'over {MAXDA:.0f} Da ({mw:.1f})')
    if m is None: reasons.append('CCD SMILES not parsed')
    row = dict(ligand_id=comp, name=cc['name'][:70], pdb_entries=';'.join(e['entry'] for e in ents), resolutions=';'.join(e['resolution_A'] for e in ents), formula=formula, mw=round(mw, 2), heavy_atoms=heavy, net_charge=Chem.GetFormalCharge(m) if m else '',
               smiles_ccd=smi, tanimoto_morgan2_to_dabogratinib=round(DataStructs.TanimotoSimilarity(dab_fp, AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048)), 3) if m else '', excluded='; '.join(reasons), retrieved=json.load(open(f'{C}/raw/ccd_{comp}.json.meta'))['retrieved'])
    if not reasons:
        mh = Chem.AddHs(m); p = AllChem.ETKDGv3(); p.randomSeed = 42; cid = AllChem.EmbedMolecule(mh, p); how = 'ETKDGv3 seed 42'
        if cid != 0: p.useRandomCoords = True; cid = AllChem.EmbedMolecule(mh, p); how = 'ETKDGv3 seed 42 + useRandomCoords (first embedding failed)'
        if cid != 0: row['build'] = 'EMBED FAILED'; rows.append(row); log.p(f'  {comp}: embedding failed'); continue
        mm = AllChem.MMFFOptimizeMolecule(mh); sdf = f'{C}/ligands/{comp}.sdf'; pq = f'{C}/ligands/{comp}.pdbqt'
        with Chem.SDWriter(sdf) as w: w.write(mh)
        cp = subprocess.run([mk, '-i', sdf, '-o', pq], capture_output=True, text=True)
        if cp.returncode != 0 or not os.path.exists(pq): row['build'] = 'MEEKO FAILED: ' + cp.stderr[-200:]; rows.append(row); log.p(f'  {comp}: Meeko failed'); continue
        row['build'] = f'{how}; MMFF return {mm}'; row['pdbqt'] = os.path.relpath(pq, L.REV); row['pdbqt_md5'] = hashlib.md5(open(pq, 'rb').read()).hexdigest(); row['rotatable_bonds_torsdof'] = sum(1 for l in open(pq) if l.startswith('TORSDOF')) and [l.split()[1] for l in open(pq) if l.startswith('TORSDOF')][0]
    rows.append(row)
keys = list(dict.fromkeys(k for r in rows for k in r))
with open(f'{C}/null_ligands.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
inc = [r for r in rows if not r['excluded'] and r.get('pdbqt')]
log.p(f'ligands: {len(rows)} distinct IDs; excluded {sum(1 for r in rows if r["excluded"])}; built {len(inc)}; failed builds {[r["ligand_id"] for r in rows if not r["excluded"] and not r.get("pdbqt")]}')
for r in rows:
    if r['excluded']: log.p(f"  EXCLUDED {r['ligand_id']} ({r['pdb_entries']}): {r['excluded']}")
log.p('  included:', ', '.join(f"{r['ligand_id']}[{r['pdb_entries']}]" for r in inc))
log.p('  Tanimoto to dabogratinib among included: ' + ', '.join(f"{r['ligand_id']} {r['tanimoto_morgan2_to_dabogratinib']}" for r in sorted(inc, key=lambda x: -x['tanimoto_morgan2_to_dabogratinib'])[:10]))
