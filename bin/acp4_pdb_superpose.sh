#!/bin/bash

set -u
#set -x #DEBUG

if [ "$#" -eq 0 ]; then
    echo "usage: acp4_pdb_superpose.sh REF.pdb MOV.pdb"
    echo "       will try to optimally superpose MOV onto REF"
    echo "       in the ph4 feature space"
    exit 1
fi

REF_PDB=$1
MOV_PDB=$2

DIR=`dirname ${MOV_PDB}`
BASE=`basename ${MOV_PDB} .pdb`
# output file: superposed PDB
SUP_PDB=${DIR}/${BASE}_sup.pdb

echo "WARNING: heuristic method" >> /dev/stderr

TMP_OUT_DIR=`mktemp -d`
#echo "DEBUG: "${TMP_OUT_DIR}
trap 'rm -rf ${TMP_OUT_DIR}' EXIT

# get ref and curr in ph4 format
obabel ${REF_PDB} -O ${TMP_OUT_DIR}/ref.sdf 2>&1 | tee ${TMP_OUT_DIR}/obabel.log
obabel ${MOV_PDB} -O ${TMP_OUT_DIR}/mov.sdf 2>&1 | tee ${TMP_OUT_DIR}/obabel.log

acp4_ph4.py -i ${TMP_OUT_DIR}/ref.sdf -o ${TMP_OUT_DIR}/ref.ph4 2>&1 | tee ${TMP_OUT_DIR}/acp4_ph4.log
acp4_ph4.py -i ${TMP_OUT_DIR}/mov.sdf -o ${TMP_OUT_DIR}/mov.ph4 2>&1 | tee ${TMP_OUT_DIR}/acp4_ph4.log

# try to find optimal transform (requires the nlopt library)
export LD_LIBRARY_PATH=/apps/nlopt/2.7.1/lib64:${LD_LIBRARY_PATH}
acp4_superpose -ref ${TMP_OUT_DIR}/ref.ph4 -cand ${TMP_OUT_DIR}/mov.ph4 \
               -oph4 ${TMP_OUT_DIR}/sup.ph4 -o ${TMP_OUT_DIR}/params.txt \
               2>&1 | tee ${TMP_OUT_DIR}/acp4_superpose.log

echo "INFO: optimal parameters:" >> /dev/stderr
tail -6 ${TMP_OUT_DIR}/params.txt >> /dev/stderr

echo "INFO: optimally superposed pdb in: "${SUP_PDB} >> /dev/stderr
acp4_pdb_move.py -i ${MOV_PDB} -ip ${TMP_OUT_DIR}/params.txt -o ${SUP_PDB}

## visualize
# pymol ${REF_PDB} ${SUP_PDB}
