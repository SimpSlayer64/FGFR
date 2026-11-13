"""validate.py - check that a docked pose lands in the FGFR3 crystal ligand site.
Takes the docked pose file as an argument. 9CD7 contains three copies of the inhibitor (A1AV2); the centroid of the pose is compared with each copy and the nearest one is reported.
Only meaningful for an FGFR3 (9CD7-derived) docking run. A centroid distance ignores orientation and is a quick check only.
Usage: python validate.py docked_pose.pdbqt"""
import sys
import numpy as np


def crystal_copies(cif='9CD7.cif', resn='A1AV2'):
    """Return {auth_chain: Nx3 array} for each copy of the ligand."""
    copies = {}
    in_atom = False
    with open(cif) as f:
        for line in f:
            if line.startswith('_atom_site.'):
                in_atom = True
            if in_atom and line.startswith('HETATM') and resn in line:
                p = line.split()
                try:
                    x, y, z = float(p[10]), float(p[11]), float(p[12])
                    chain = p[-5]  # auth_asym_id
                    copies.setdefault(chain, []).append([x, y, z])
                except (ValueError, IndexError):
                    pass
    return {k: np.array(v) for k, v in copies.items()}


def docked_first_pose(path):
    coords = []
    with open(path) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            if line.startswith('ENDMDL'):
                break
    return np.array(coords)


def main():
    if len(sys.argv) != 2:
        raise SystemExit('usage: python validate.py <docked_pose.pdbqt>')
    docked_path = sys.argv[1]

    copies = crystal_copies()
    print(f'Crystal ligand copies found: {len(copies)} '
          f'({", ".join(f"{k}:{len(v)}at" for k, v in copies.items())})')

    docked = docked_first_pose(docked_path)
    print(f'Docked pose atoms: {len(docked)}  (from {docked_path})')
    docked_center = docked.mean(axis=0)

    dists = {k: float(np.linalg.norm(v.mean(axis=0) - docked_center))
             for k, v in copies.items()}
    nearest = min(dists, key=dists.get)
    d = dists[nearest]

    print(f'Docked centroid : {docked_center.round(2)}')
    for k, dv in sorted(dists.items(), key=lambda kv: kv[1]):
        print(f'  vs crystal copy {k} centroid {copies[k].mean(0).round(2)} : {dv:.2f} A')
    print(f'Nearest crystal copy: {nearest}  ({d:.2f} A)')

    if d < 5.0:
        print('VALIDATION PASSED (pose overlaps a crystal TYRA-300 copy)')
    elif d < 10.0:
        print('PARTIAL: pose near a crystal copy but not well overlapped')
    else:
        print('VALIDATION FAILED: pose far from every crystal TYRA-300 copy '
              '(is this actually an FGFR3/9CD7 docking run?)')


if __name__ == '__main__':
    main()
