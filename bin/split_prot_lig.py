#!/usr/bin/env python3

# split XXXX.pdb into XXXX_lig.pdb and XXXX_rec.pdb
# several ligands might end-up in XXXX_lig.pdb

import os, sys
from prody import *

pdb_fn = sys.argv[1]
base_fn = pdb_fn.removesuffix('.pdb')

try:
    pdb = parsePDB(pdb_fn)
except:
    print("prody: cannot parse: %s" % pdb_fn, file=sys.stderr)
    sys.exit(1)
else:
    ligand = pdb.select('(not protein) and (not water) and (not ion)')
    if ligand:
        writePDB(base_fn + '_lig.pdb', ligand)
    receptor = pdb.select('protein')
    if receptor:
        writePDB(base_fn + '_rec.pdb', receptor)
