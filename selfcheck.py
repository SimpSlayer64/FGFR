"""selfcheck.py - runs the main steps in a scratch folder and compares with the stored values (the source of each expected value is noted next to it).
Needs vina, meeko, rdkit, gemmi and Open Babel. Set FGFR3_DOCK to a folder that holds 9CD7.cif, 9CD5.cif, receptor_fgfr3_v2.pdbqt and receptor_fgfr1_v2.pdbqt
(the supplied receptors are used for docking because Open Babel's polar-hydrogen placement is not bitwise reproducible). Writes only to a temporary folder.
Usage: FGFR3_DOCK=/path/to/data python selfcheck.py"""
import os, sys, json, csv, shutil, subprocess, tempfile
B = os.path.dirname(os.path.abspath(__file__)); D = os.environ['FGFR3_DOCK']; T = tempfile.mkdtemp(prefix='fgfr3_selfcheck_'); fails = []
def run(*a, **k): r = subprocess.run([sys.executable, *a], cwd=T, capture_output=True, text=True, env={**os.environ, 'FGFR3_DOCK': T, 'PYTHONDONTWRITEBYTECODE': '1'}); r.check_returncode() if not k.get('soft') else None; return r
def check(name, ok, detail): print(('PASS ' if ok else 'FAIL ') + name + ': ' + detail); (None if ok else fails.append(name))
atoms = lambda f: [l for l in open(f) if l.startswith(('ATOM', 'HETATM'))]
for f in ('9CD7.cif', '9CD5.cif'): shutil.copy(f'{D}/{f}', T)
# 1 receptors (box centres exact; atom counts as in the committed files; hydrogens may differ in position)
run(f'{B}/build_receptors.py'); bc = json.load(open(f'{T}/box_centers.json'))
check('box centres', bc == {'fgfr3': [26.76, -7.51, 4.11], 'fgfr1': [-19.2, 2.13, 8.12]}, str(bc))
for t in ('fgfr3', 'fgfr1'): check(f'{t} receptor atoms', len(atoms(f'{T}/receptor_{t}_v2.pdbqt')) == len(atoms(f'{D}/receptor_{t}_v2.pdbqt')), f'{len(atoms(f"{T}/receptor_{t}_v2.pdbqt"))} rebuilt vs {len(atoms(f"{D}/receptor_{t}_v2.pdbqt"))} committed')
# 2 ligands (deterministic: the rebuilt dabogratinib must equal the committed file atom for atom)
run(f'{B}/build_ligands.py'); check('dabogratinib ligand', atoms(f'{T}/dabogratinib_corr.pdbqt') == atoms(f'{D}/dabogratinib_corr.pdbqt'), 'ATOM lines identical to the committed dabogratinib_corr.pdbqt')
# 3 docking, seed 11, exhaustiveness 8, box 30, committed receptors (expected: replicate runs, replicate_scores.csv, seed 11: FGFR3 -9.418, FGFR1 -8.966)
for t, exp in (('fgfr3', -9.418), ('fgfr1', -8.966)):
    shutil.copy(f'{D}/receptor_{t}_v2.pdbqt', f'{T}/committed_{t}.pdbqt')
    r = run(f'{B}/dock.py', '--ligand', 'dabogratinib_corr.pdbqt', '--receptor', f'committed_{t}.pdbqt', '--center', *map(str, bc[t]), '--seed', '11', '--out', f'pose_{t}.pdbqt')
    got = float([l for l in r.stdout.splitlines() if l.startswith('Best binding')][0].split(':')[1].split()[0])
    check(f'dock {t} seed 11', abs(got - exp) <= 0.005, f'{got} vs {exp}')
# 4 RMSD (expected: rmsd.csv, stage 2f_fiveseed, seed 11, model 1, crystal copy A: FGFR3 1.352)
sys.path.insert(0, B); import build_ligands as BL
run(f'{B}/p2_rmsd.py', '--cif', '9CD7.cif', '--chain', 'A', '--ligand-resn', 'A1AV2', '--template-smiles', BL.COMPOUNDS['dabogratinib'][0], '--poses', 'pose_fgfr3.pdbqt', '--out', 'rmsd.csv')
v = [float(x['rmsd']) for x in csv.DictReader(open(f'{T}/rmsd.csv')) if x['model'] == '1' and x['crystal_copy'] == 'A'][0]
check('RMSD FGFR3 model 1, copy A', abs(v - 1.352) <= 0.001, f'{v} vs 1.352')
shutil.rmtree(T, ignore_errors=True); print('RESULT:', 'all checks passed' if not fails else f'FAILED: {fails}'); sys.exit(1 if fails else 0)
