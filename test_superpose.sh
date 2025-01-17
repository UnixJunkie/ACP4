#!/bin/bash

# move ligand pdb
./bin/acp4_pdb_move.py -i data/3rfm_xtal_lig.pdb -o moved_lig.pdb -r 1.5,1.0,0.5 -t 1.,5.,10.

# get ref and curr in ph4 format
