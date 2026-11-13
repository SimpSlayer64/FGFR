"""Prepares every selected receptor (no docking): protein residues of one chain only (no ligand, water, ions), `obabel <chain>.pdb -O <out>.pdbqt -xr -p 7.4`,
loads it into Vina (set_receptor only) and computes the hinge superposition, the box centre and the superposed 9CD7 ligand centroid for each receptor.
Hydrogens already present in the deposited file are removed first; alternative conformations are removed with gemmi remove_alternative_conformations (first conformer kept);
the chain is renamed A in the written PDB. Hinge = backbone N, CA, C, O of FGFR3-numbered residues 556-559 (FGFR1 562-565) = 16 atoms.
Reference for the box: 9CD7 chain A (FGFR3 receptors) or 9CD5 chain A (FGFR1 receptors); box centres from box_centers.json.
Usage: PYTHONDONTWRITEBYTECODE=1 python p3b_prep.py"""
import sys, os, csv, io, json, subprocess, contextlib, importlib.util, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_lib as L
import numpy as np, gemmi
with contextlib.redirect_stdout(io.StringIO()):
    spec = importlib.util.spec_from_file_location('sel', f'{L.REV}/phase3b/p3b_select.py'); S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
log = L.Log('p3b_prep'); OUT = f'{L.REV}/phase3b/ensemble'; RD = f'{OUT}/receptors'
ob = subprocess.run(['obabel', '-V'], capture_output=True, text=True).stdout.strip(); log.p('obabel:', ob)
BB = S.BB; HING = S.HINGE; BC = json.load(open(f'{L.DOCK}/box_centers.json'))
def kabsch(P, Q):      # maps Q onto P
    pc, qc = P.mean(0), Q.mean(0); H = np.matmul((Q - qc).T, P - pc); U, Sv, Vt = np.linalg.svd(H); d = np.sign(np.linalg.det(np.matmul(Vt.T, U.T)))
    R = np.matmul(np.matmul(Vt.T, np.diag([1, 1, d])), U.T); return R, pc - np.matmul(R, qc)
def hinge(m): return np.array([[a.pos.tolist() for a in m[h] if a.name == n][0] for h in HING for n in BB])
def fit(rec_h, ref_h):
    R, t = kabsch(rec_h, ref_h); mapped = np.matmul(ref_h, R.T) + t; return R, t, float(np.sqrt(((mapped - rec_h) ** 2).sum(1).mean()))
# references: hinge coordinates and the 9CD7 chain A ligand centroid
REF = {}
for ent in ('9CD7', '9CD5'):
    p = f'{L.DOCK}/{ent}.cif'; st = gemmi.read_structure(p); st.setup_entities(); b = gemmi.cif.read(p).sole_block()
    iso, m, ch, _ = S.chain_map(st, b, 'A'); REF[ent] = dict(h=hinge(m), lig=np.array([a.pos.tolist() for r in ch if r.name == 'A1AV2' for a in r]))
C7 = REF['9CD7']['lig'].mean(0); log.p('9CD7 chain A A1AV2 centroid (heavy atoms)', np.round(C7, 3), '| box_centers.json fgfr3', BC['fgfr3'], '| 9CD5 chain A centroid', np.round(REF['9CD5']['lig'].mean(0), 3), '| box_centers.json fgfr1', BC['fgfr1'])
rows = []
def process(entry, chain, iso, res, mut, kind, path_cif, ready_pdbqt=None):
    d = dict(entry=entry, chain=chain, isoform=iso, resolution_A=res, mutations=mut, kind=kind)
    try:
        st = gemmi.read_structure(path_cif); st.setup_entities(); b = gemmi.cif.read(path_cif).sole_block()
        if kind == 'ensemble': st.remove_alternative_conformations(); st.remove_hydrogens()
        _, m, ch, _ = S.chain_map(st, b, chain); h = hinge(m)
        R7, t7, r7 = fit(h, REF['9CD7']['h']); d['hinge_rmsd_to_9CD7'] = round(r7, 3)
        d['ref9CD7_ligand_centroid_in_receptor_frame'] = ' '.join(f'{x:.3f}' for x in (np.matmul(C7, R7.T) + t7))
        refname = '9CD7' if iso == 'FGFR3' else '9CD5'
        Rb, tb, rb = fit(h, REF[refname]['h']); d['box_reference'] = refname; d['hinge_rmsd_to_box_reference'] = round(rb, 3)
        c0 = np.array(BC['fgfr3' if iso == 'FGFR3' else 'fgfr1']); c = np.matmul(c0, Rb.T) + tb; d['cx'], d['cy'], d['cz'] = [round(float(x), 3) for x in c]
        if kind == 'reference':
            c = c0; d['cx'], d['cy'], d['cz'] = [float(x) for x in c0]; d['box_note'] = 'box_centers.json as docked in the standard runs (reference, no transformation)'
            if entry == '9CD7': d['ref9CD7_ligand_centroid_in_receptor_frame'] = ' '.join(f'{x:.3f}' for x in C7)
        d['box_centre_minus_9CD7_ligand_centroid_A'] = round(float(np.linalg.norm(c - np.array([float(x) for x in d['ref9CD7_ligand_centroid_in_receptor_frame'].split()]))), 3)
        if ready_pdbqt: pq = ready_pdbqt
        else:
            prot = gemmi.Structure(); prot.add_model(gemmi.Model('1')); nch = gemmi.Chain('A')
            for res_ in ch:
                if S.aa(res_): nch.add_residue(res_)
            prot[0].add_chain(nch); pdb = f'{RD}/{entry}_{chain}_chainA.pdb'; prot.write_pdb(pdb); pq = f'{RD}/{entry}_{chain}_v2.pdbqt'
            cp = subprocess.run(['obabel', pdb, '-O', pq, '-xr', '-p', '7.4'], capture_output=True, text=True); open(f'{RD}/{entry}_{chain}_obabel.log', 'w').write(cp.stdout + cp.stderr)
            d['obabel_returncode'] = cp.returncode; d['n_residues_written'] = len(nch)
        n = sum(1 for l in open(pq) if l.startswith('ATOM')); hd = sum(1 for l in open(pq) if l.startswith('ATOM') and l[77:79].strip() == 'HD')
        d['receptor_file'] = os.path.relpath(pq, L.DOCK); d['receptor_md5'] = hashlib.md5(open(pq, 'rb').read()).hexdigest(); d['n_atoms'] = n; d['n_HD'] = hd
        from vina import Vina
        Vina(sf_name='vina', cpu=1, verbosity=0).set_receptor(pq); d['vina_load_ok'] = True
    except Exception as ex: d['vina_load_ok'] = False; d['error'] = f'{type(ex).__name__}: {ex}'
    rows.append(d)
    log.p(f"{entry} {chain} {iso} {kind}: ok={d.get('vina_load_ok')} atoms {d.get('n_atoms')} HD {d.get('n_HD')} hinge rmsd -> 9CD7 {d.get('hinge_rmsd_to_9CD7')} -> {d.get('box_reference')} {d.get('hinge_rmsd_to_box_reference')} box ({d.get('cx')}, {d.get('cy')}, {d.get('cz')}) {d.get('error', '')}")
process('9CD7', 'A', 'FGFR3', '', '', 'reference', f'{L.DOCK}/9CD7.cif', f'{L.DOCK}/receptor_fgfr3_v2.pdbqt')
process('9CD5', 'A', 'FGFR1', '', '', 'reference', f'{L.DOCK}/9CD5.cif', f'{L.DOCK}/receptor_fgfr1_v2.pdbqt')
for r in csv.DictReader(open(f'{OUT}/../receptor_list.csv')):
    process(r['entry'], r['chain'], r['isoform'], r['resolution_A'], r['mutations'], 'ensemble', f"{S.P3}/raw/pdb/{r['entry']}.cif")
keys = list(dict.fromkeys(k for d in rows for k in d))
with open(f'{OUT}/receptors.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
bad = [d for d in rows if not d.get('vina_load_ok')]
log.p(f'\nreceptors {len(rows)} (2 reference + {len(rows) - 2} ensemble); Vina load failures: {[(d["entry"], d.get("error")) for d in bad]}')
for k in ('hinge_rmsd_to_9CD7', 'hinge_rmsd_to_box_reference', 'box_centre_minus_9CD7_ligand_centroid_A'):
    v = [d[k] for d in rows if k in d and d['kind'] == 'ensemble']; log.p(f'{k}: ensemble min {min(v)} median {float(np.median(v)):.3f} max {max(v)}')
