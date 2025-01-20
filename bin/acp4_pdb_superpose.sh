#!/bin/bash

set -u
set -x #DEBUG

# FBR: TODO usage message in case no args
#      - test
#      - check all scripts are installed by opam for acp4
#      - install latest ACP4

REF_PDB=$1
MOV_PDB=$2

DIR=`dirname ${MOV_PDB}`
BASE=`basename ${MOV_PDB}`
# output file: superposed PDB
SUP_PDB=${DIR}/${BASE}_sup.pdb

echo "WARNING: heuristic method" >> /dev/stderr

TMP_OUT_DIR=`mktemp -d`
# echo ${TMP_OUT_DIR} # DEBUG
trap 'rm -rf ${TMP_OUT_DIR}' EXIT

# get ref and curr in ph4 format
obabel ${REF_PDB} -O ${TMP_OUT_DIR}/ref.sdf 2>&1 | tee ${TMP_OUT_DIR}/obabel.log
obabel ${MOV_PDB} -O ${TMP_OUT_DIR}/mov.sdf 2>&1 | tee ${TMP_OUT_DIR}/obabel.log

acp4_ph4.py -i ${TMP_OUT_DIR}/ref.sdf -o ${TMP_OUT_DIR}/ref.ph4 2>&1 | tee ${TMP_OUT_DIR}/acp4_ph4.log
acp4_ph4.py -i ${TMP_OUT_DIR}/mov.sdf -o ${TMP_OUT_DIR}/mov.ph4 2>&1 | tee ${TMP_OUT_DIR}/acp4_ph4.log

# try to find optimal transform
# (acp4_superpose requires NLopt library)
export LD_LIBRARY_PATH=/apps/nlopt/2.7.1/lib64:${LD_LIBRARY_PATH}
acp4_superpose -ref ${TMP_OUT_DIR}/ref.ph4 -cand ${TMP_OUT_DIR}/mov.ph4 \
               -o ${TMP_OUT_DIR}/sup.ph4 2>&1 | tee ${TMP_OUT_DIR}/acp4_superpose.log

# apply it
# WARNING: prefix each float w/ a space
# so that negative numbers are not considered CLI options
# FBR: get transform params from last line in prev. command
acp4_pdb_move.py -i ${MOV_PDB} -o opt.pdb \
                 -r " -1.525810, -0.559073, 0.977297" \
                 -to " 7.737021, -33.434615, -32.893111"

# # visualize if DEBUG
# pymol data/3rfm_xtal_lig.pdb opt.pdb
