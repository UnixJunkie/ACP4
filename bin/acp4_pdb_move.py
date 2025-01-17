#!/usr/bin/env python3

# rotate then translate a PDB
# rotation given by angles in radians: rx, ry, rz
# translation relative to the current position: tx, ty, tz

import argparse
import numpy as np
import sys
import typing

from prody import moveAtoms, parsePDB, writePDB
from prody.measure import calcCenter
from prody.measure.transform import Transformation
from scipy.spatial.transform import Rotation as Rot

origin = np.array([0., 0., 0.])
no_trans = np.array([0., 0., 0.])

def rot_from_rx_ry_rz(rx: float, ry: float, rz: float):
    rx = Rot.from_euler('x', rx).as_matrix()
    ry = Rot.from_euler('y', ry).as_matrix()
    rz = Rot.from_euler('z', rz).as_matrix()
    r = np.matmul(np.matmul(rz, ry), rx)
    return Transformation(r, no_trans)

def read_three_floats(line: str) -> tuple[float]:
    x, y, z = line.strip().split(',')
    return (float(x), float(y), float(z))

if __name__ == "__main__":
    # CLI options parsing -----------------------------------------------------
    parser = argparse.ArgumentParser(
        description = "rotate then translate a pdb file")
    parser.add_argument("-i", metavar = "input.pdb", dest = "input_fn",
                        help = "input file")
    parser.add_argument("-o", metavar = "output.pdb", dest = "output_fn",
                        help = "output file")
    parser.add_argument("-r", metavar = "rx,ry,rz", dest = "r_xyz",
                        help = "angular rotations", default='0.,0.,0.')
    parser.add_argument("-t", metavar = "tx,ty,tz", dest = "t_xyz",
                        help = "relative translations", default='0.,0.,0.')
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    # -------------------------------------------------------------------------
    input_fn: str = args.input_fn
    output_fn: str = args.output_fn
    rx, ry, rz = read_three_floats(args.r_xyz)
    tx, ty, tz = read_three_floats(args.t_xyz)
    # 1) center PDB
    atoms = parsePDB(input_fn)
    previous_center = calcCenter(atoms)
    print('prev_center: %s' % previous_center)
    new_center = previous_center + np.array([tx, ty, tz])
    print('next_center: %s' % new_center)
    moveAtoms(atoms, to=origin)
    # 2) rotate around Ox, Oy then Oz
    rot = rot_from_rx_ry_rz(rx, ry, rz)
    rot.apply(atoms)
    # 3) translate to final position
    moveAtoms(atoms, to=new_center)
    # 4) write out pdb
    writePDB(output_fn, atoms)
