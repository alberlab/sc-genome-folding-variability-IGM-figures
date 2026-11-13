#!/bin/bash
#$ -N extfig2c
#$ -cwd
#$ -o /u/project/falber/wyeast/H1_hESC_five_stages_Dekker/extfig2_occ06/logs/extfig2c.$JOB_ID.$TASK_ID.out
#$ -j y
#$ -l h_data=6G,h_rt=2:00:00
#$ -t 1-4
set -euo pipefail
unset LD_LIBRARY_PATH
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
source /u/project/xjzhou/wyeast/miniconda3/etc/profile.d/conda.sh
conda activate igm_schic_exp_env
R=/u/project/falber/wyeast/H1_hESC_five_stages_Dekker
D=$R/extfig2_occ06
read CELL HSS A B C < <(sed -n "${SGE_TASK_ID}p" $R/fig5_occ06/manifest.tsv)
echo "task $SGE_TASK_ID cell $CELL start $(date) on $(hostname)"
LAB=$R/locale_4cell_occ06_N565/locale_occ06_N565/locale_labels_${CELL}.npz
SPK=$R/g3d_sf_4models_occ06/$CELL/speckle_clusters
python $D/extfig2c_locale_structures.py "$CELL" "$HSS" "$A" "$B" "$C" $LAB \
    $R/fig5_occ06/inputs/${CELL}_occ06_locale_N565.tsv $SPK $D/out/$CELL
cd $D/out/$CELL && tar czf ${CELL}_extfig2c_locale_cmm_pdb.tar.gz cmm_pdb_locale speckles ${CELL}_extfig2c.cxc ${CELL}_locale_bead_counts_per_structure.tsv
echo "EXTFIG2C_DONE $CELL $(date)"
