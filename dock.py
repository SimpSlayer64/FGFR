"""dock.py - single Vina docking run with explicit arguments and an output file that is never overwritten silently.

Example:
  python dock.py --ligand dabogratinib_corr.pdbqt --receptor receptor_fgfr3_v2.pdbqt --center 26.76 -7.51 4.11 --seed 11 --out pose_fgfr3.pdbqt"""
import argparse
import os
from vina import Vina


def main():
    ap = argparse.ArgumentParser(description='Single Vina docking run.')
    ap.add_argument('--ligand',   required=True, help='prepared ligand .pdbqt')
    ap.add_argument('--receptor', required=True, help='prepared receptor .pdbqt')
    ap.add_argument('--center',   required=True, type=float, nargs=3,
                    metavar=('X', 'Y', 'Z'), help='box center in Angstroms')
    ap.add_argument('--box',      type=float, nargs=3, default=[30, 30, 30],
                    metavar=('X', 'Y', 'Z'), help='box size (default 30 30 30)')
    ap.add_argument('--out',      required=True,
                    help='output poses .pdbqt (must be unique per run)')
    ap.add_argument('--exhaustiveness', type=int, default=8)
    ap.add_argument('--n-poses',        type=int, default=5)
    ap.add_argument('--seed',           type=int, default=42,
                    help='fixed RNG seed so runs are reproducible (default 42)')
    ap.add_argument('--overwrite', action='store_true',
                    help='allow overwriting --out if it already exists')
    args = ap.parse_args()

    if os.path.exists(args.out) and not args.overwrite:
        raise SystemExit(
            f'refusing to overwrite existing {args.out} (pass --overwrite to force)')

    v = Vina(sf_name='vina', seed=args.seed, verbosity=1)
    v.set_receptor(args.receptor)
    v.set_ligand_from_file(args.ligand)
    v.compute_vina_maps(center=list(args.center), box_size=list(args.box))
    v.dock(exhaustiveness=args.exhaustiveness, n_poses=args.n_poses)
    v.write_poses(args.out, n_poses=args.n_poses, overwrite=True)

    energies = v.energies()
    print('Docking complete')
    print(f'Ligand   : {args.ligand}')
    print(f'Receptor : {args.receptor}')
    print(f'Center   : {args.center}   Box: {args.box}   Seed: {args.seed}')
    print(f'Output   : {args.out}')
    print(f'Best binding affinity: {energies[0][0]:.3f} kcal/mol')
    print(f'All poses: {[round(e[0], 3) for e in energies]}')


if __name__ == '__main__':
    main()
