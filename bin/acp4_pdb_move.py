#!/usr/bin/env python3

import prody, sys

if __name__ == "__main__":

input_fn = sys.argv[1]
rx = sys.argv[2]
ry = sys.argv[3]
rz = sys.argv[4]
tx = sys.argv[5]
ty = sys.argv[6]
tz = sys.argv[7]

# 1) center PDB
# 2) rotate around Ox, Oy then Oz
# 3) translate to final position
