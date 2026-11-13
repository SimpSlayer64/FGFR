"""Writes docking manifests (no docking). Usage: python p2_manifests.py anchor|m2f|e2|c2"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_lib as L
import csv
which = sys.argv[1]
bc = L.box_centres()
ORDER = [r['compound'] for r in csv.DictReader(open(f'{L.REV}/compounds.csv'))]
V2 = {'fgfr3': 'receptor_fgfr3_v2.pdbqt', 'fgfr1': 'receptor_fgfr1_v2.pdbqt'}
def lig(n): return f'work/ligands/{n}_rebuilt.pdbqt'
def j(group, ligand, receptor, centre, seed, tag='', cpu=1):
    # jobs on the corrected receptors with the canonical ligand at the standard conditions always live in group 'replicates'
    shared = ligand.startswith('work/ligands/') and ligand.endswith('_rebuilt.pdbqt') and 'oldisomer' not in ligand and receptor in V2.values() and seed in L.SEEDS5
    return L.make_job('replicates' if shared else group, ligand, receptor, centre, box=30, exh=8, seed=seed, cpu=cpu, conf='conf0', tag=tag)
jobs = []
if which == 'anchor':
    for n in ORDER:
        for t in ('fgfr3', 'fgfr1'): jobs.append(L.make_job('anchor', lig(n), V2[t], bc[t], seed=42, cpu=1, conf='conf0', tag='anchor seed42'))
elif which == 'm2f':
    for s in L.SEEDS5:
        for t in ('fgfr3', 'fgfr1'): jobs.append(j('replicates', lig('dabogratinib'), V2[t], bc[t], s, '2f/2c'))
elif which == 'e2':
    old3, old1 = ('receptor_clean.pdbqt', [26.76, -7.51, 4.11]), ('receptor_fgfr1.pdbqt', [-34.20, 29.56, 38.04])
    arms = {'a_oldfile': 'ligand_9CD7.pdbqt', 'b_oldrebuilt': 'work/ligands/dabogratinib_oldisomer_rebuilt.pdbqt', 'c_correct': lig('dabogratinib')}
    for s in L.SEEDS5:
        for arm, lg in arms.items():
            jobs.append(j('attrib', lg, V2['fgfr3'], bc['fgfr3'], s, f'2e {arm} FGFR3 v2'))
            jobs.append(j('attrib', lg, V2['fgfr1'], bc['fgfr1'], s, f'2e {arm} FGFR1 9CD5'))
            jobs.append(j('attrib', lg, old1[0], old1[1], s, f'2e {arm} FGFR1 old3GQI'))
            jobs.append(j('attrib', lg, old3[0], old3[1], s, f'2e-ext {arm} FGFR3 old receptor_clean'))
        cl = lig('custom1')
        jobs.append(j('attrib', cl, V2['fgfr3'], bc['fgfr3'], s, '2e custom1 FGFR3 v2'))
        jobs.append(j('attrib', cl, V2['fgfr1'], bc['fgfr1'], s, '2e custom1 FGFR1 9CD5'))
        jobs.append(j('attrib', cl, old1[0], old1[1], s, '2e custom1 FGFR1 old3GQI'))
        jobs.append(j('attrib', cl, old3[0], old3[1], s, '2e-ext custom1 FGFR3 old receptor_clean'))
elif which == 'c2':
    for n in ORDER:
        for s in L.SEEDS5:
            for t in ('fgfr3', 'fgfr1'): jobs.append(j('replicates', lig(n), V2[t], bc[t], s, '2c'))
seen = {}; 
for x in jobs: seen[x['id']] = x
jobs = list(seen.values())
os.makedirs(f'{L.REV}/manifests', exist_ok=True)
L.write_manifest(f'{L.REV}/manifests/{which}.csv', jobs)
new = [x for x in jobs if not os.path.exists(L.meta_path(x))]
print(f'{which}: {len(jobs)} unique jobs, {len(new)} not yet run -> manifests/{which}.csv')
