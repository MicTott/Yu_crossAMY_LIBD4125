_#!/bin/bash
#$ -cwd
#$ -l mem_free=32G,h_vmem=32G,h_fsize=300G
#$ -N Yu_crossAMY_LIBD4125_cellranger_aggr
#$ -o logs/cellranger_aggr.$TASK_ID.txt
#$ -e logs/cellranger_aggr.$TASK_ID.txt
#$ -m e
#$ -t 1-3
#$ -tc 3

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## load CellRanger
module load cellranger/7.0.0

## List current modules for reproducibility
module list

## Locate file
SAMPLE=$(awk "NR==${SGE_TASK_ID}" 04_cellranger_aggr.txt)
echo "Processing sample ${SAMPLE}"
date

## Run CellRanger aggr
cellranger aggr --id=${SAMPLE}_aggr --csv=./aggr_files/${SAMPLE}.csv

## Move output
echo "Moving data to new location"
date
mv ${SAMPLE}_aggr /dcs04/lieber/marmaypag/Yu_crossAMY_LIBD4125/processed-data/03_cellranger/

echo "**** Job ends ****"
date
