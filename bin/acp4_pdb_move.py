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

def read_last_six_params(fn):
    lines = open(fn).readlines()
    lines = list(map(lambda s: s.strip(), lines))
    n = len(lines)
    if n == 6:
        return (float(lines[0]),
                float(lines[1]),
                float(lines[2]),
                float(lines[3]),
                float(lines[4]),
                float(lines[5]))
    elif n == 12:
        return (float(lines[6]),
                float(lines[7]),
                float(lines[8]),
                float(lines[9]),
                float(lines[10]),
                float(lines[11]))
    else:
        print('acp4_pdb_move.py: read_last_six_params: file %s has %d lines' %
              (fn, n), file=sys.stderr)
        exit(1)

def read_center_from_file(center_fn):
    line = open(center_fn).readline()
    line = line.strip()
    toks = line.split(',')
    x = float(toks[0])
    y = float(toks[1])
    z = float(toks[2])
    return [x, y, z]

if __name__ == "__main__":
    # CLI options parsing -----------------------------------------------------
    parser = argparse.ArgumentParser(
        description = "rotate then translate a pdb file")
    parser.add_argument("-i", metavar = "input.pdb", dest = "input_fn",
                        help = "input PDB to move")
    parser.add_argument("-ip", metavar = "input.txt", dest = "in_params_fn",
                        help = "optimal parameters input file")
    parser.add_argument("-o", metavar = "output.pdb", dest = "output_fn",
                        help = "moved PDB output")
    parser.add_argument('--center', dest='center_fn', type=str,
                        help = "binding-site x,y,z coordinates in a file",
                        default = None)
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    # -------------------------------------------------------------------------
    input_fn: str = args.input_fn
    in_params_fn: str = args.in_params_fn
    output_fn: str = args.output_fn
    rx, ry, rz, dx, dy, dz = read_last_six_params(in_params_fn)
    center_fn = args.center_fn
    bs_center = None
    if center_fn != None:
        bs_center = read_center_from_file(center_fn)
    # 1) center PDB
    atoms = parsePDB(input_fn)
    current_center = calcCenter(atoms)
    cx, cy, cz = current_center[0], current_center[1], current_center[2]
    # output current center in case we need to apply the same transform later
    # to another pdb
    print('prev_center: %f %f %f' % (cx, cy, cz))
    new_center = np.array([dx, dy, dz])
    print('next_center: %s' % new_center, file=sys.stderr)
    moveAtoms(atoms, to=origin)
    if bs_center != None:
        print('moving BS center to origin', file=sys.stderr)
        bs_displacement = current_center - bs_center
        print('bs_displacement: %s' % bs_displacement, file=sys.stderr)
        moveAtoms(atoms, by=bs_displacement)
    # 2) rotate around Ox, Oy then Oz
    rot = rot_from_rx_ry_rz(rx, ry, rz)
    rot.apply(atoms)
    # 3) translate to final position
    moveAtoms(atoms, by=new_center)
    # 4) write out pdb
    writePDB(output_fn, atoms)
