"""Selects the receptor ensemble from the table of PDB entries (phase3/pdb_fgfr1_fgfr3_entries.csv). No docking.
Criteria:
  C1 resolution <= 2.5 A; 9CD7 and 9CD5 are left out.
  C2 a drug-sized ligand (MW > 250, not an additive, not a nucleotide: ACP ANP ADP ATP AGS AMP) in the same chain with a heavy atom within 4.0 A of the backbone N or O of the hinge alanine (FGFR3 A558 = FGFR1 A564).
  C3 no missing residues within 6 A of the ligand: the pocket set is every residue with a heavy atom within 6.0 A of the inhibitor copies in any chain of 9CD7 or 9CD5 (FGFR1 number - 6 = FGFR3 number);
     every position of the set must be observed in the chain with N, CA, C and O.
  C0 the chain is an FGFR1 or FGFR3 kinase domain: its construct sequence aligns to the 9CD5 (FGFR1) or 9CD7 (FGFR3) construct with >= 70 percent identity over >= 200 aligned residues.
     Residues are numbered through this alignment, not through _struct_ref_seq.
  One receptor per entry: the lowest chain letter that meets C0, C2 and C3.
Usage: PYTHONDONTWRITEBYTECODE=1 python p3b_select.py"""
import os, sys, csv, json, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_lib as L
import numpy as np, gemmi
log = L.Log('p3b_select')
OUT = f'{L.REV}/phase3b'; P3 = f'{L.REV}/phase3'
FGFR_UNI = {'P11362': 'FGFR1', 'P22607': 'FGFR3'}; OFF = {'FGFR1': 6, 'FGFR3': 0}
NUC = {'ACP', 'ANP', 'ADP', 'ATP', 'AGS', 'AMP'}
HINGE = [556, 557, 558, 559]; HINGE_ALA = 558; BB = ('N', 'CA', 'C', 'O')
RES_MAX, HINGE_CUT, POCKET_CUT = 2.5, 4.0, 6.0
REFS = {'9CD7': f'{L.DOCK}/9CD7.cif', '9CD5': f'{L.DOCK}/9CD5.cif'}

def aa(res):
    i = gemmi.find_tabulated_residue(res.name); return bool(i and i.is_amino_acid())
def nw(s, t, m=3, mm=-1, gp=-4):
    """global alignment, free end gaps; returns list of (i, j) aligned index pairs """
    n, k = len(s), len(t); S = np.zeros((n + 1, k + 1)); T = np.zeros((n + 1, k + 1), dtype=int)
    for i in range(1, n + 1): T[i, 0] = 1
    for j in range(1, k + 1): T[0, j] = 2
    for i in range(1, n + 1):
        for j in range(1, k + 1):
            c = [S[i-1, j-1] + (m if s[i-1] == t[j-1] else mm), S[i-1, j] + (0 if j == k else gp), S[i, j-1] + (0 if i == n else gp)]
            a = int(np.argmax(c)); S[i, j] = c[a]; T[i, j] = a
    i, j, p = n, k, []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and T[i, j] == 0: p.append((i-1, j-1)); i -= 1; j -= 1
        elif i > 0 and (j == 0 or T[i, j] == 1): i -= 1
        else: j -= 1
    return p[::-1]
def entity_seq(st, b, chain):
    """one-letter construct sequence of the chain's entity (full SEQRES, not only observed residues)"""
    for ent in st.entities:
        if chain in [sub for sub in ent.subchains] or True:
            pass
    ch = [c for c in st[0] if c.name == chain][0]
    sc = ch[0].subchain if len(ch) else None
    for ent in st.entities:
        if ent.entity_type == gemmi.EntityType.Polymer and sc in ent.subchains:
            return ent, ''.join(gemmi.find_tabulated_residue(n).one_letter_code.upper() if gemmi.find_tabulated_residue(n) and gemmi.find_tabulated_residue(n).one_letter_code != ' ' else 'X' for n in ent.full_sequence)
    return None, ''
REFSEQ = {}
for _ent, _path in REFS.items():
    _st = gemmi.read_structure(_path); _st.setup_entities(); _b = _st.make_mmcif_document().sole_block() if False else gemmi.cif.read(_path).sole_block()
    _e, _seq = entity_seq(_st, _b, 'A')
    _t = _b.find('_struct_ref_seq.', ['pdbx_strand_id', 'db_align_beg', 'seq_align_beg']); _row = [r for r in _t if r[0] == 'A'][0]
    _start = int(_row[1]) - (int(_row[2]) - 1)
    REFSEQ[_ent] = (_seq, _start - (6 if _ent == '9CD5' else 0), 'FGFR3' if _ent == '9CD7' else 'FGFR1')   # start = FGFR3 number of construct residue 1
ALN_CACHE = {}
def chain_map(st, b, chain):
    """residue -> FGFR3 numbering through sequence alignment to 9CD7 / 9CD5. Returns isoform, map, chain, identity table."""
    ch = [c for c in st[0] if c.name == chain][0]
    ent, seq = entity_seq(st, b, chain)
    if not seq: return None, {}, ch, {}
    best = None; info = {}
    for ref, (rseq, rstart, riso) in REFSEQ.items():
        key = (seq, ref)
        if key not in ALN_CACHE: ALN_CACHE[key] = nw(seq, rseq)
        al = ALN_CACHE[key]; ident = sum(1 for i, j in al if seq[i] == rseq[j]); info[riso] = (round(ident / max(1, len(al)), 3), len(al))
        if best is None or ident > best[0]: best = (ident, ref, al, rstart, riso)
    ident, ref, al, rstart, riso = best
    if not (info[riso][0] >= 0.70 and info[riso][1] >= 200): return None, {}, ch, info
    idx2fn = {i + 1: rstart + j for i, j in al}                      # label_seq (1-based) -> FGFR3 number
    m = {}
    for res in ch:
        if res.label_seq is None or not aa(res): continue
        if res.label_seq in idx2fn: m[idx2fn[res.label_seq]] = res
    return riso, m, ch, info
def has_bb(res): return all(any(a.name == n for a in res) for n in BB)
def xyz(res, names=None): return np.array([a.pos.tolist() for a in res if names is None or a.name in names])
def ligands(ch, comps):
    return [r for r in ch if r.name in comps and r.name not in NUC and not aa(r)]

# ---- reference pocket set (union over all TYRA-300 copies of 9CD7 and 9CD5)
pocket = {}; ref_info = {}
for ent, path in REFS.items():
    st = gemmi.read_structure(path); st.setup_entities(); b = gemmi.cif.read(path).sole_block()
    for c in st[0]:
        iso, m, ch, _inf = chain_map(st, b, c.name)
        for lig in [r for r in ch if r.name == 'A1AV2']:
            L_xyz = xyz(lig); near = []
            for fn, res in m.items():
                P = xyz(res)
                if np.linalg.norm(P[:, None, :] - L_xyz[None, :, :], axis=2).min() <= POCKET_CUT: near.append(fn)
            for fn in near: pocket.setdefault(fn, set()).add(f'{ent}:{c.name}')
            ref_info[(ent, c.name)] = (iso, len(near), sorted(near))
            log.p(f'reference {ent} chain {c.name} ({iso}): {len(near)} residues within {POCKET_CUT} A of A1AV2; hinge {[ (h, m[h].name if h in m else "unobs") for h in HINGE]}')
POCKET = sorted(pocket)
log.p(f'pocket set (union over {len(ref_info)} TYRA-300 copies): {len(POCKET)} positions (FGFR3 numbering): {POCKET}')
# residue-name identity check across isoforms (FGFR1 number - 6)
names = {}
for ent, path in REFS.items():
    st = gemmi.read_structure(path); st.setup_entities(); b = gemmi.cif.read(path).sole_block()
    iso, m, ch, _inf = chain_map(st, b, 'A'); names[ent] = {k: v.name for k, v in m.items()}
log.p('positions of the pocket set with different residue names in 9CD7 A and 9CD5 A (FGFR3 numbering):', [(p, names['9CD7'].get(p), names['9CD5'].get(p)) for p in POCKET if names['9CD7'].get(p) != names['9CD5'].get(p)])

# ---- candidates: every chain with a listed ligand, for every entry of the table
rows = list(csv.DictReader(open(f'{P3}/pdb_fgfr1_fgfr3_entries.csv')))
meta = {}
for r in rows: meta.setdefault(r['entry'], r)
out = []
for e, r in meta.items():
    comps = set(r['ligand_comp_ids'].split(';'))
    base = dict(entry=e, isoform_phase3_list=r['isoform'], resolution_A=r['resolution_A'], ligand_comp_ids=r['ligand_comp_ids'], mutations=r['mutations_in_entity'], title=r['title'])
    base['C1_resolution_le_2.5'] = bool(r['resolution_A']) and float(r['resolution_A']) <= RES_MAX
    base['excluded_reference'] = e in REFS
    if not base['C1_resolution_le_2.5'] or base['excluded_reference']:
        out.append(dict(base, chain='', passes_all=False, note='not examined further (C1 or reference)')); continue
    path = f'{P3}/raw/pdb/{e}.cif'
    st = gemmi.read_structure(path); st.setup_entities(); b = gemmi.cif.read(path).sole_block()
    chains = [c for c in st[0] if any(x.name in comps and not aa(x) for x in c)]
    if not chains: out.append(dict(base, chain='', passes_all=False, note='no chain holds a ligand of the candidate table')); continue
    for c in chains:
        d = dict(base, chain=c.name, C0_fgfr_kinase=False, C2_ligand_in_ATP_site=False, ligand_resname='', ligand_heavy_atoms='', hinge_Ala_min_dist_A='', C3_pocket_complete=False,
                 pocket_missing_positions='', n_pocket_positions=len(POCKET), note='')
        iso, m, ch, info = chain_map(st, b, c.name)
        d['identity_FGFR1_FGFR3'] = f"{info.get('FGFR1', '')} / {info.get('FGFR3', '')}"
        d['isoform'] = iso or ''
        d['C0_fgfr_kinase'] = iso is not None
        if iso is None: d['note'] = 'not an FGFR1/FGFR3 kinase domain by sequence identity; '
        else:
            best = None
            if HINGE_ALA in m and has_bb(m[HINGE_ALA]):
                hb = xyz(m[HINGE_ALA], ('N', 'O'))
                for lig in ligands(ch, comps):
                    dd = float(np.linalg.norm(xyz(lig)[:, None, :] - hb[None, :, :], axis=2).min())
                    if best is None or dd < best[0]: best = (dd, lig)
            else: d['note'] += 'hinge Ala558-equivalent not observed; '
            if best:
                d['hinge_Ala_min_dist_A'] = round(best[0], 2); d['ligand_resname'] = best[1].name; d['ligand_heavy_atoms'] = sum(1 for a in best[1] if a.element.name != 'H')
                d['C2_ligand_in_ATP_site'] = best[0] <= HINGE_CUT
            elif not d['note']: d['note'] += 'no non-nucleotide ligand of the candidate table in this chain; '
            miss = [p for p in POCKET if p not in m or not has_bb(m[p])]
            d['pocket_missing_positions'] = ','.join(map(str, miss)); d['C3_pocket_complete'] = not miss
            mism = [(p, m[p].name, names['9CD7' if iso == 'FGFR3' else '9CD5'].get(p)) for p in POCKET if p in m and m[p].name != names['9CD7' if iso == 'FGFR3' else '9CD5'].get(p)]
            d['pocket_residue_differences_vs_reference'] = ';'.join(f'{p}:{a}vs{b_}' for p, a, b_ in mism)
            d['hinge_names'] = ''.join(gemmi.find_tabulated_residue(m[h].name).one_letter_code.upper() if h in m else '-' for h in HINGE)
        d['passes_all'] = bool(d['C0_fgfr_kinase'] and d['C2_ligand_in_ATP_site'] and d['C3_pocket_complete'])
        out.append(d)
keys = list(dict.fromkeys(k for d in out for k in d))
with open(f'{OUT}/receptor_candidates_all_chains.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(out)
sel = {}
for d in sorted(out, key=lambda x: (x['entry'], x['chain'])):
    if d['passes_all'] and d['entry'] not in sel: sel[d['entry']] = d
S = sorted(sel.values(), key=lambda x: (x['isoform'], x['entry']))
with open(f'{OUT}/receptor_list.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['entry', 'chain', 'isoform', 'resolution_A', 'ligand_resname', 'ligand_heavy_atoms', 'hinge_Ala_min_dist_A', 'hinge_names', 'pocket_residue_differences_vs_reference', 'mutations', 'title']); w.writeheader()
    for d in S: w.writerow({k: d.get(k, '') for k in w.fieldnames})
ents = collections.defaultdict(list)
for d in out: ents[d['entry']].append(d)
log.p(f'\ncandidate table: {len(rows)} chain rows, {len(ents)} entries')
c1 = [e for e, v in ents.items() if v[0]['C1_resolution_le_2.5']]
log.p(f'C1 resolution <= 2.5: {len(c1)} entries; reference entries among them (left out): {[e for e in c1 if e in REFS]}')
c1n = [e for e in c1 if e not in REFS]
nofgfr = [e for e in c1n if not any(d.get('C0_fgfr_kinase') for d in ents[e])]
log.p(f'C0 not an FGFR1/FGFR3 kinase by sequence: {nofgfr}: ' + str({e: ents[e][0].get('identity_FGFR1_FGFR3') for e in nofgfr}))
c2 = [e for e in c1n if any(d.get('C0_fgfr_kinase') and d.get('C2_ligand_in_ATP_site') for d in ents[e])]
c3 = [e for e in c2 if any(d['passes_all'] for d in ents[e])]
log.p(f'C2 drug-like ligand in ATP site: {len(c2)} entries (fail C2, among FGFR kinases: {sorted(set(c1n) - set(nofgfr) - set(c2))})')
log.p(f'C3 no missing pocket residues (in a chain that also meets C2): {len(c3)} entries; fail C3 only: {sorted(set(c2) - set(c3))}')
log.p(f'selected: {len(S)} receptors: FGFR1 {sum(1 for d in S if d["isoform"] == "FGFR1")}, FGFR3 {sum(1 for d in S if d["isoform"] == "FGFR3")}')
for d in S: log.p(f"  {d['entry']} chain {d['chain']} {d['isoform']} {d['resolution_A']} A {d['ligand_resname']} ({d['ligand_heavy_atoms']} heavy) hinge {d['hinge_names']} dist {d['hinge_Ala_min_dist_A']} {d['mutations'] or 'WT'} | pocket differences: {d['pocket_residue_differences_vs_reference'] or 'none'}")
for e in sorted(set(c2) - set(c3)): log.p('  C3 fail', e, [(d['chain'], d['pocket_missing_positions']) for d in ents[e] if d.get('C2_ligand_in_ATP_site')])
for e in sorted(set(c1n) - set(nofgfr) - set(c2)): log.p('  C2 fail', e, [(d['chain'], d['ligand_resname'], d['hinge_Ala_min_dist_A'], d['note']) for d in ents[e]])
