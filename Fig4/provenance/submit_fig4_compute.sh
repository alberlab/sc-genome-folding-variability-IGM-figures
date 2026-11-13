#!/bin/bash
#$ -N fig4cmp
#$ -cwd
#$ -o /u/project/falber/wyeast/H1_hESC_five_stages_Dekker/fig4_occ06/logs/fig4cmp.$JOB_ID.$TASK_ID.out
#$ -j y
#$ -l h_data=8G,h_rt=4:00:00
#$ -t 1-4
set -euo pipefail
unset LD_LIBRARY_PATH
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
source /u/project/xjzhou/wyeast/miniconda3/etc/profile.d/conda.sh
conda activate igm_schic_exp_env
D=/u/project/falber/wyeast/H1_hESC_five_stages_Dekker/fig4_occ06
read CELL HSS A B C < <(sed -n "${SGE_TASK_ID}p" $D/manifest.tsv)
echo "task $SGE_TASK_ID cell $CELL start $(date) on $(hostname)"
[[ -s $HSS ]] || { echo "missing $HSS"; exit 1; }
python $D/fig4_compute.py "$CELL" "$HSS" "$A" "$B" "$C" /u/scratch/w/wyeast/ct_strictcis_2p5r/msk_strict/$CELL/top30/chunks $D/out/$CELL
cd $D/out/$CELL && tar czf ${CELL}_cmm_pdb_chrom.tar.gz cmm_pdb_chrom ${CELL}_nucleus_ellipsoid.cxc && tar czf ${CELL}_cmm_pdb_shannon.tar.gz cmm_pdb_shannon ${CELL}_nucleus_ellipsoid.cxc
echo "FIG4_DONE $CELL $(date)"
