#!/usr/bin/env python3

# extract PDB lines corresponding to a user-provided pocket definition
# you must use a conda environment w/ prody
# a pocket is specified via a text file only made of
# ^CHAIN:char<SPACE>RESNUM:int$ lines

import os, sys
from prody import *

def load_pocket_definition(pocket_fn):
    res = []
    for line in open(pocket_fn).readlines():
        strip = line.strip()
        split = line.split()
        chain = str(split[0])
        if len(chain) != 1:
            print("FATAL: get_pocket.py: load_pocket_definition: strange chain: %s" % line,
                  file=sys.stderr)
            exit(1)
        resnum = int(split[1])
        if resnum <= 0:
            print("FATAL: get_pocket.py: load_pocket_definition: strange resnum: %s" % line,
                  file=sys.stderr)
            exit(1)
        chain_res = (chain, resnum)
        res.append(chain_res)
    return res

def group_by_chain(residues):
    res = {}
    for chain, resnum in residues:
        prev = res.get(chain, [])
        prev.append(resnum)
        res[chain] = prev
    return res

# FBR: use proper CLI handling maybe one day
#      or, at least usage message in case not enough params
pdb_fn = sys.argv[1]
pocket_fn = sys.argv[2]
output_fn = sys.argv[3]

try:
    pdb = parsePDB(pdb_fn)
except:
    print("prody: cannot parse: %s" % pdb_fn, file=sys.stderr)
    sys.exit(1)
else:
    pocket_residues = load_pocket_definition(pocket_fn)
    chain_residues = group_by_chain(pocket_residues)
    with open(output_fn, 'w') as output:
        for chain, residues in chain_residues.items():
            resnums = list(map(str, residues))
            resnums_str = ' '.join(resnums)
            selection = 'protein and chain %c and resnum %s' % (chain, resnums_str)
            atoms = pdb.select(selection)
            if atoms:
                writePDBStream(output, atoms)
