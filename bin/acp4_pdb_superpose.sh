#!/bin/bash

set -u
#set -x #DEBUG

if [ "$#" -eq 0 ]; then
    echo "usage: acp4_pdb_superpose.sh REF.pdb MOV.pdb WHO.pdb"
    echo "       will try to optimally superpose MOV onto REF"
    echo "       in the ph4 feature space; then apply the same"
    echo "       transformation to WHO"
    echo "       REF: binding-site 1"
    echo "       MOV: binding-site 2"
    echo "       WHO: whole PDB for MOV"
    exit 1
fi

REF_PDB=$1
MOV_PDB=$2
ALL_PDB=$3

DIR=`dirname ${MOV_PDB}`
BASE=`basename ${MOV_PDB} .pdb`
# output file: superposed binding-site
SUP_PDB=${DIR}/${BASE}_sup_BS.pdb
WHO_PDB=${DIR}/${BASE}_sup_whole.pdb

echo "WARNING: heuristic method" >> /dev/stderr

TMP_OUT_DIR=`mktemp -d`
#echo "DEBUG: "${TMP_OUT_DIR}
trap 'rm -rf ${TMP_OUT_DIR}' EXIT

# get ref and curr in ph4 format
obabel ${REF_PDB} -O ${TMP_OUT_DIR}/ref.sdf 2>&1 >> ${TMP_OUT_DIR}/obabel.log
obabel ${MOV_PDB} -O ${TMP_OUT_DIR}/mov.sdf 2>&1 >> ${TMP_OUT_DIR}/obabel.log

# conda setup
source /apps/miniconda3/4.10.3/etc/profile.d/conda.sh
conda activate /apps/conda_envs/rdkit_prody

acp4_ph4.py --permissive -i ${TMP_OUT_DIR}/ref.sdf -o ${TMP_OUT_DIR}/ref.ph4 \
            2>&1 >  ${TMP_OUT_DIR}/acp4_ph4.log
acp4_ph4.py --permissive -i ${TMP_OUT_DIR}/mov.sdf -o ${TMP_OUT_DIR}/mov.ph4 \
            2>&1 >> ${TMP_OUT_DIR}/acp4_ph4.log

# try to find optimal transform (requires the nlopt library)
export LD_LIBRARY_PATH=/apps/nlopt/2.7.1/lib64:${LD_LIBRARY_PATH}
acp4_superpose -ref ${TMP_OUT_DIR}/ref.ph4 -cand ${TMP_OUT_DIR}/mov.ph4 \
               -oph4 ${TMP_OUT_DIR}/sup.ph4 -o ${TMP_OUT_DIR}/params.txt \
               2>&1 | tee ${TMP_OUT_DIR}/acp4_superpose.log

echo "INFO: optimal parameters:" >> /dev/stderr
tail -6 ${TMP_OUT_DIR}/params.txt >> /dev/stderr

echo "INFO: optimally superposed binding-site in: "${SUP_PDB} >> /dev/stderr
acp4_pdb_move.py -i ${MOV_PDB} -ip ${TMP_OUT_DIR}/params.txt -o ${SUP_PDB} \
                 2>&1 > ${TMP_OUT_DIR}/acp4_pdb_move.log
BS_CENTER=`grep -F 'prev_center:' ${TMP_OUT_DIR}/acp4_pdb_move.log | awk '{print $2","$3","$4}'`
echo 'BS_CENTER: '${BS_CENTER} >> /dev/stderr
echo ${BS_CENTER} > ${TMP_OUT_DIR}/center.txt

# make the whole PDB undergo the optimal transformation
echo "INFO: optimally superposed whole pdb in: "${WHO_PDB} >> /dev/stderr
acp4_pdb_move.py -i ${ALL_PDB} -ip ${TMP_OUT_DIR}/params.txt \
                 -o ${WHO_PDB} --center ${TMP_OUT_DIR}/center.txt \
                 2>&1 >> ${TMP_OUT_DIR}/acp4_pdb_move.log

## visualize
# pymol ${REF_PDB} ${SUP_PDB}
