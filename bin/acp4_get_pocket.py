#!/usr/bin/env python3

# extract PDB lines corresponding to a user-provided pocket definition
# you must use a conda environment w/ prody
# a pocket is specified via a text file only made of
# ^CHAIN:char<SPACE>RESNUM:int$ lines

import os, sys
from prody import parsePDB, writePDBStream

def load_pocket_definition(pocket_fn):
    res = []
    for line in open(pocket_fn).readlines():
        strip = line.strip()
        split = line.split()
        chain = str(split[0])
        if len(chain) != 1:
            print("FATAL: acp4_get_pocket.py: load_pocket_definition: strange chain: %s" % line,
                  file=sys.stderr)
            exit(1)
        resnum = int(split[1])
        if resnum <= 0:
            print("FATAL: acp4_get_pocket.py: load_pocket_definition: strange resnum: %s" % line,
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
output_fn = sys.argv[3] # this one is either a .pdb file (default) or a .sdf
                        # in which case obabel will be called to do pdb2sdf

pdb_output_fn = ''
sdf_output_fn = ''

convert_pdb_to_sdf = False

if output_fn.endswith('.sdf'):
    convert_pdb_to_sdf = True
    pdb_output_fn = output_fn.replace('.sdf', '.pdb')
    sdf_output_fn = output_fn

if output_fn.endswith('.pdb'):
    convert_pdb_to_sdf = False
    pdb_output_fn = output_fn
    sdf_output_fn = ''

try:
    pdb = parsePDB(pdb_fn)
except:
    print("prody: cannot parse: %s" % pdb_fn, file=sys.stderr)
    sys.exit(1)
else:
    pocket_residues = load_pocket_definition(pocket_fn)
    chain_residues = group_by_chain(pocket_residues)
    with open(pdb_output_fn, 'w') as output:
        for chain, residues in chain_residues.items():
            resnums = list(map(str, residues))
            resnums_str = ' '.join(resnums)
            selection = 'protein and chain %c and resnum %s' % (chain, resnums_str)
            atoms = pdb.select(selection)
            if atoms:
                writePDBStream(output, atoms)
    if sdf_output_fn != '':
        # try calling obabel for pdb2sdf conversion
        cmd = 'obabel %s -O %s' % (pdb_output_fn, sdf_output_fn)
        print("INFO: acp4_get_pocket.py: running: %s" % cmd,
              file=sys.stderr)
        os.system(cmd)
