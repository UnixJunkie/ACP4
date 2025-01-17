#!/bin/bash

# test pdb translations

# # tx
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.,0.,0. -t 1.,0.,0. -o tx1.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.,0.,0. -t 2.,0.,0. -o tx2.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.,0.,0. -t 3.,0.,0. -o tx3.pdb

# pymol test.pdb tx1.pdb tx2.pdb tx3.pdb

# # ty
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.,0.,0. -t 0.,1.,0. -o ty1.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.,0.,0. -t 0.,2.,0. -o ty2.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.,0.,0. -t 0.,3.,0. -o ty3.pdb

# pymol test.pdb ty1.pdb ty2.pdb ty3.pdb

# # tz
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.,0.,0. -t 0.,0.,1. -o tz1.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.,0.,0. -t 0.,0.,2. -o tz2.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.,0.,0. -t 0.,0.,3. -o tz3.pdb

# pymol test.pdb tz1.pdb tz2.pdb tz3.pdb

# # test rotations
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.2,0.,0. -t 1.,0.,0. -o tx1.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.4,0.,0. -t 2.,0.,0. -o tx2.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.6,0.,0. -t 3.,0.,0. -o tx3.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.8,0.,0. -t 4.,0.,0. -o tx4.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 1.0,0.,0. -t 5.,0.,0. -o tx5.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 1.2,0.,0. -t 6.,0.,0. -o tx6.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 1.4,0.,0. -t 7.,0.,0. -o tx7.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 1.6,0.,0. -t 8.,0.,0. -o tx8.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 1.8,0.,0. -t 9.,0.,0. -o tx9.pdb

# pymol test.pdb tx*.pdb

# # test rotations
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.,0.2,0. -t 0.,1,0. -o tx1.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.,0.4,0. -t 0.,2,0. -o tx2.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.,0.6,0. -t 0.,3,0. -o tx3.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.,0.8,0. -t 0.,4,0. -o tx4.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.,1.0,0. -t 0.,5,0. -o tx5.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.,1.2,0. -t 0.,6,0. -o tx6.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.,1.4,0. -t 0.,7,0. -o tx7.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.,1.6,0. -t 0.,8,0. -o tx8.pdb
# ./bin/acp4_pdb_move.py -i test.pdb -r 0.,1.8,0. -t 0.,9,0. -o tx9.pdb

# pymol test.pdb tx*.pdb
