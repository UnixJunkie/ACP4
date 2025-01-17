#!/usr/bin/env python3

# rotate then translate a PDB
# rotation given by angles in radians: rx, ry, rz
# translation relative to the current position: tx, ty, tz

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

if __name__ == "__main__":
    # FBR: proper CLI
    input_fn: str = sys.argv[1]
    rx = float(sys.argv[2])
    ry = float(sys.argv[3])
    rz = float(sys.argv[4])
    tx = float(sys.argv[5])
    ty = float(sys.argv[6])
    tz = float(sys.argv[7])
    output_fn = sys.argv[8]
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
