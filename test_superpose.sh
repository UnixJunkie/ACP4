#!/bin/bash

# move ligand pdb
./bin/acp4_pdb_move.py -i data/3rfm_xtal_lig.pdb -o moved_lig.pdb -r " 1.5, 1.0, 0.5" -t " 1.0, 5.0, 10.0"

# get ref and curr in ph4 format
obabel data/3rfm_xtal_lig.pdb -O ref.sdf
obabel moved_lig.pdb          -O moved.sdf

# project in ph4 space
./bin/acp4_ph4.py -i ref.sdf   -o ref.ph4
./bin/acp4_ph4.py -i moved.sdf -o moved.ph4

# try to find optimal transform
./superpose -ref ref.ph4 -cand moved.ph4 -o sup.ph4

# apply it
# WARNING: prefix each float w/ a space
# so that negative numbers are not considered CLI options
./bin/acp4_pdb_move.py -i moved_lig.pdb -o opt.pdb \
                       -r " -1.525810, -0.559073, 0.977297" \
                       -to " 7.737021, -33.434615, -32.893111"

# visualize
pymol data/3rfm_xtal_lig.pdb opt.pdb

# get binding-site for CFF in PDB:3rfm
acp4_get_pocket.py data/3rfm.pdb data/3rfm_pocket.txt data/3rfm_pocket.pdb
# get binding-site for CFF in PDB:5mzp
acp4_get_pocket.py data/5mzp.pdb data/5mzp_pocket.txt data/5mzp_pocket.pdb

./bin/acp4_pdb_superpose.sh data/3rfm_pocket.pdb data/5mzp_pocket.pdb data/5mzp.pdb

# visualize
pymol data/3rfm.pdb data/3rfm_pocket.pdb data/5mzp.pdb data/5mzp_pocket.pdb \
      data/5mzp_pocket_sup_BS.pdb data/5mzp_pocket_sup_whole.pdb
