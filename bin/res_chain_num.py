#!/usr/bin/env python3

# extract (chain, res_name) pairs from a PDB file for HETATM lines
# only print a res_name once (i.e. only 1st time given ligand is seen)

import os, sys
from prody import *

pdb_fn = sys.argv[1]
pdb = pdb_fn.removesuffix(".pdb")

# PDB line format
# 0                17  2123
#"HETATM    1  N   VWW A 210      ..."
#res_name: off:17 len:3
#chain: off:21 len:1
#res_num: off:22 len:4
def triplet_of_line(line):
    res_name = line[17:21].strip()
    chain = line[21:22]
    res_num = line[22:26].strip()
    return (res_name, chain, res_num)

s = set()
for line in open(pdb_fn, 'r').readlines():
    if line.startswith('HETATM'):
        name, chain, _num = triplet_of_line(line)
        chain_lig = (chain, name)
        s.add(chain_lig)

seen = set()
for chain, res_name in s:
    if res_name not in seen:
        print('%s_%s_%s' % (pdb, chain, res_name))
        seen.add(res_name)
