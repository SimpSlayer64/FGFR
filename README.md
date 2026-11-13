FGFR3 / FGFR1 docking scripts

Scripts that prepare the receptors and ligands, dock with Vina, measure pose accuracy and run the receptor ensemble for the FGFR3 versus FGFR1 selectivity index (SI = best FGFR3 score minus best FGFR1 score, kcal/mol).

SOFTWARE VERSIONS
- Python 3.11.15
- Vina 1.2.7 (Python API)
- Meeko 0.7.1
- Open Babel 3.2.0
- RDKit 2026.03.1
- PyMOL 3.1.8 (figures only; no script here needs it)
- gemmi 0.7.5, numpy 1.26.4

INPUT DATA
- data/receptor_fgfr3_v2.pdbqt: FGFR3 receptor (PDB 9CD7, chain A).
- data/receptor_fgfr1_v2.pdbqt: FGFR1 receptor (PDB 9CD5, chain A).
- data/box_centers.json: box centres (30 x 30 x 30 A boxes).
- data/dabogratinib_corr.pdbqt: dabogratinib ligand (also produced by build_ligands.py).
- 9CD7.cif and 9CD5.cif are not included; download them from the Protein Data Bank.
Open Babel places polar hydrogens in a way that is not bitwise reproducible. build_receptors.py therefore gives receptors that are equivalent to, but not identical with, the supplied files. Use the supplied receptors to reproduce the stored scores.

FOLDERS
- FGFR3_DOCK (environment variable, default: current folder): holds 9CD7.cif, 9CD5.cif, the receptor files, box_centers.json and the ligand files.
- FGFR3_REV (default: $FGFR3_DOCK/work): working folder for runs/, logs/ and manifests/. The ensemble scripts expect phase3/pdb_fgfr1_fgfr3_entries.csv (table of PDB entries) and phase3/raw/pdb/<ID>.cif (mmCIF files) here and write phase3b/.
- Keep all scripts in one folder; they import p2_lib.py from their own folder. Set PYTHONDONTWRITEBYTECODE=1 to avoid bytecode folders.

RECEPTOR COMMAND (build_receptors.py, run for each chain-A file, protein residues only)
  obabel chain.pdb -O receptor.pdbqt -xr -p 7.4

DOCKING SETTINGS
Vina scoring function, box 30 x 30 x 30 A, exhaustiveness 8, 5 output poses, score = best pose. Seeds: 11, 22, 33, 44, 55 for the main runs; 42 for the single runs and for every conformer embedding of the base ligands; 101 and 202 for the two extra conformers of the sensitivity test. Thread count does not change the scores; the seed does.

RUN ORDER
1. build_receptors.py: receptor_fgfr3_v2.pdbqt, receptor_fgfr1_v2.pdbqt, box_centers.json (needs 9CD7.cif, 9CD5.cif in the current folder).
2. build_ligands.py: the 21 ligands as <name>_corr.sdf and <name>_corr.pdbqt (RDKit ETKDGv3 seed 42, MMFF94, Meeko).
3. phase1j_rebuild_ligands.py <docking_dir> <work_dir>: rebuilds the ligands, compares them with step 2 and writes the alternative-conformer set.
4. p2_manifests.py anchor|m2f|e2|c2: writes the docking manifests to manifests/.
5. p2_run_jobs.py manifests/<name>.csv: runs a manifest, one Vina docking per job, resumable; writes runs/<group>/ (pose .pdbqt and .json per job) and docking_runs.csv.
6. p2_rmsd_stage.py 2f|2b_seed42|2b_replicates (uses p2_rmsd.py): RMSD of poses against the crystal ligand copies; writes rmsd.csv. p2_rmsd_sanity.py checks the RMSD method; p2_rmsd_summary.py summarises rmsd.csv.
7. p2_replicates.py: five-seed summary (mean SI, SD, 95 percent t-interval with t = 2.776); writes replicate_scores.csv and replicate_summary.csv.
8. p2_sensitivity.py select, then run the manifest written to manifests/ with p2_run_jobs.py, then p2_sensitivity.py summary: sensitivity.csv.
9. p2g_prep.py: receptor re-preparation test (work/prep/, manifests/g2.csv).
10. validate.py <pose.pdbqt>: centroid check of one FGFR3 pose against the 9CD7 crystal ligand copies.
11. p3b_select.py: receptor_list.csv (ensemble selection).
12. p3b_prep.py: ensemble receptors and box centres (phase3b/ensemble/receptors/, receptors.csv).
13. p3b_dock.py: docks dabogratinib into the reference pair and the ensemble (seeds 11, 22, 33).
14. p3b_analyze.py: ensemble_summary.csv, ensemble_pairs.csv, ensemble_receptor_scores.csv, ensemble_per_run.csv, ENSEMBLE_RESULT.txt.
15. p3c_null_build.py: builds the set of unrelated ATP-site ligands for the same docking.
Single docking: dock.py --ligand L.pdbqt --receptor R.pdbqt --center X Y Z --seed S --out pose.pdbqt (prints the best score to 3 decimals).

EXPECTED OUTPUT (check with selfcheck.py)
python selfcheck.py with FGFR3_DOCK pointing to a folder that holds 9CD7.cif, 9CD5.cif and the four files of data/ runs steps 1, 2, a seed 11 docking of dabogratinib and the RMSD, in a temporary folder. Expected:
- box_centers.json: fgfr3 [26.76, -7.51, 4.11], fgfr1 [-19.2, 2.13, 8.12].
- receptor atom counts 2652 (FGFR3) and 2744 (FGFR1).
- dabogratinib_corr.pdbqt ATOM lines identical to data/dabogratinib_corr.pdbqt.
- seed 11 best score: -9.418 (FGFR3) and -8.966 (FGFR1), within 0.005.
- RMSD of the FGFR3 top pose to crystal copy A: 1.352 A.
The ensemble scripts need the PDB table and mmCIF files and are not covered by selfcheck.py.

VERIFIED
These checks were run from a copy of this folder placed in an empty folder outside the project (selfcheck.py and the ligand test), with the software versions listed above.
- Seed 11 docking of dabogratinib (data/dabogratinib_corr.pdbqt, 30 A box, exhaustiveness 8, supplied receptors), best score in kcal/mol, against the stored value:
    FGFR3 (9CD7): -9.418 (stored -9.418)
    FGFR1 (9CD5): -8.966 (stored -8.966)
- Redock RMSD of the FGFR3 top pose (seed 11) to crystal ligand copy A, symmetry-aware: 1.352 A (stored 1.352 A).
- build_receptors.py: box centres identical to data/box_centers.json; atom counts 2652 (FGFR3) and 2744 (FGFR1), equal to the supplied receptor files. Hydrogen positions differ slightly from the supplied files, as stated above.
- build_ligands.py: all 21 compounds were built in an empty folder. The InChIKey of each output SDF equals the InChIKey listed for that compound in the supplementary table of isomeric SMILES and InChIKeys (21 of 21 equal, 0 mismatches). The dabogratinib ATOM lines equal those of data/dabogratinib_corr.pdbqt.
- Every script compiles, and p2_lib.py imports and runs its statistics function.

NOT RUN
The scripts below were not run in the clean copy. Each needs inputs that are not shipped here.
- phase1j_rebuild_ligands.py: compounds.csv (table of the 21 compounds with SMILES and in_paper_20 column) and the existing <name>_corr.sdf / <name>_corr.pdbqt files of the docking folder.
- p2_manifests.py: compounds.csv and the rebuilt ligand files in work/ligands/ (output of the step above).
- p2_run_jobs.py: a manifest from p2_manifests.py or p2_sensitivity.py, and all ligand and receptor files named in it.
- p2_rmsd_stage.py: compounds.csv, 9CD7.cif and 9CD5.cif, the pose files in work/runs/ (output of p2_run_jobs.py) and, for stages 2f and 2b_seed42, the earlier single-run pose files docked_<name>_<receptor>_corr.pdbqt.
- p2_rmsd_summary.py: rmsd.csv (output of p2_rmsd_stage.py) and pose_drift_corr.csv (table of centroid drift of the single-run poses); it writes rmsd_summary.csv.
- p2_rmsd_sanity.py: 9CD7.cif.
- p2_replicates.py: compounds.csv, rmsd.csv, paper_values.csv (table of earlier single-run scores) and all five-seed run records in work/runs/replicates/.
- p2_sensitivity.py: replicate_summary.csv, the base run records, the extra conformer ligands it builds, and sensitivity run records.
- p2g_prep.py: 9CD7.cif and 9CD5.cif (it writes the chain-A files itself) and the ligand file work/ligands/dabogratinib_rebuilt.pdbqt.
- validate.py: a docked FGFR3 pose file and 9CD7.cif in the current folder.
- p3b_select.py: phase3/pdb_fgfr1_fgfr3_entries.csv (table of FGFR1 and FGFR3 PDB entries) and the mmCIF file of every entry in phase3/raw/pdb/.
- p3b_prep.py: receptor_list.csv (output of p3b_select.py), the same mmCIF files, data/box_centers.json and Open Babel.
- p3b_dock.py: the prepared ensemble receptors (output of p3b_prep.py), work/ligands/dabogratinib_rebuilt.pdbqt and the stored replicate run records for the reference check.
- p3b_analyze.py: the run records and pose files written by p3b_dock.py.
- p3c_null_build.py: receptor_list.csv, compounds.csv, internet access to the RCSB data API for the ligand records, and the mmCIF files.
- The ensemble scripts above (p3b_*, p3c_null_build.py) are not covered by selfcheck.py.
