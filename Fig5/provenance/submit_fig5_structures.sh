#!/bin/bash
#$ -N fig5str
#$ -cwd
#$ -o /u/project/falber/wyeast/H1_hESC_five_stages_Dekker/fig5_occ06/logs/fig5str.$JOB_ID.$TASK_ID.out
#$ -j y
#$ -l h_data=6G,h_rt=2:00:00
#$ -t 1-4
set -euo pipefail
unset LD_LIBRARY_PATH
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
source /u/project/xjzhou/wyeast/miniconda3/etc/profile.d/conda.sh
conda activate igm_schic_exp_env
D=/u/project/falber/wyeast/H1_hESC_five_stages_Dekker/fig5_occ06
read CELL HSS A B C < <(sed -n "${SGE_TASK_ID}p" $D/manifest.tsv)
echo "task $SGE_TASK_ID cell $CELL start $(date) on $(hostname)"
SPK=/u/scratch/w/wyeast/g3d_sf_4models_occ06/$CELL/speckle_clusters
python $D/fig5_structures.py "$CELL" "$HSS" "$A" "$B" "$C" $D/inputs/${CELL}_occ06_distCT_top30_2p5r_stats.tsv \
    $D/inputs/${CELL}_occ06_locale_N565.tsv $SPK $D/out/$CELL
cd $D/out/$CELL && tar czf ${CELL}_fig5_cmm_pdb.tar.gz cmm_pdb_IQR cmm_pdb_SD cmm_pdb_entropy speckles ${CELL}_fig5.cxc
# archive the per-structure speckle clusters (scratch only so far)
ARC=/u/project/falber/wyeast/H1_hESC_five_stages_Dekker/g3d_sf_4models_occ06/$CELL/speckle_clusters
mkdir -p $ARC && rsync -a $SPK/ $ARC/ && echo "archived $(ls $ARC | wc -l) speckle cluster files to $ARC"
echo "FIG5_DONE $CELL $(date)"
