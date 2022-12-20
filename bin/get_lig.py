#!/usr/bin/env python3

# extract ligand and receptor from given PDB_CHAIN_LIG triplet; e.g. 1yom_A_P01

import os, sys
from prody import *

pdb_chain_lig = sys.argv[1]

tokens = pdb_chain_lig.split('_')
pdb = tokens[0]
pdb_fn = pdb_chain_lig + ".pdb"
chain_letter = tokens[1]
lig_res = tokens[2]

# # debug
# print(chain_letter)
# print(lig_res)

# # D/L the PDB file; overwriting any existing one
# os.system('wget https://files.rcsb.org/download/%s.pdb -O %s' % \
#           (pdb, pdb_fn))
try:
    pdb = parsePDB(pdb_fn)
except:
    print("prody: cannot parse: %s" % pdb_fn, file=sys.stderr)
    sys.exit(1)
else:
    ligand = pdb.select('(chain %s) and (resname %s) and (not protein) and \
    (not water) and (not ion)' % (chain_letter, lig_res))
    if ligand:
        writePDB(pdb_chain_lig + '_lig.pdb', ligand)
    receptor = pdb.select('(chain %s) and protein' % chain_letter)
    if receptor:
        writePDB(pdb_chain_lig + '_rec.pdb', receptor)
