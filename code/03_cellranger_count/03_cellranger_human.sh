_#!/bin/bash
#$ -cwd
#$ -l mem_free=16G,h_vmem=16G,h_fsize=200G
#$ -pe local 4
#$ -N Yu_crossAMY_LIBD4125_cellranger
#$ -o logs/cellranger_human.$TASK_ID.txt
#$ -e logs/cellranger_human.$TASK_ID.txt
#$ -m e
#$ -t 1-14
#$ -tc 4

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
SAMPLE=$(awk "NR==${SGE_TASK_ID}" 03_cellranger_human.txt)
echo "Processing sample ${SAMPLE}"
date

## Run CellRanger
cellranger count --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-GRCh38-3.0.0 \
    --fastqs=/dcs04/lieber/marmaypag/Yu_crossAMY_LIBD4125/raw-data/FASTQ/${SAMPLE} \
    --sample=${SAMPLE} \
    --jobmode=local \
    --localcores=4 \
    --localmem=64

## Move output
echo "Moving data to new location"
date
mkdir -p /dcs04/lieber/marmaypag/Yu_crossAMY_LIBD4125/processed-data/03_cellranger/
mv ${SAMPLE} /dcs04/lieber/marmaypag/Yu_crossAMY_LIBD4125/processed-data/03_cellranger/

echo "**** Job ends ****"
date