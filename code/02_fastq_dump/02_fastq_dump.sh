#!/bin/bash
#$ -cwd
#$ -l mem_free=64G,h_vmem=64G,h_fsize=300G
#$ -N Yu_crossAMY_LIBD4125_FASTQ
#$ -o ./logs/FASTQ.$TASK_ID.txt
#$ -e ./logs/FASTQ.$TASK_ID.txt
#$ -m e
#$ -t 1-19
#$ -tc 10


# GEO location
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195445

# ######### Human ############
# SRR17762195
# SRR17762196
# SRR17762197
# SRR17762198
# SRR17762199
# SRR17762200
# SRR17762201
# SRR17762202
# SRR17762203
# SRR17762204
# SRR17762205
# SRR17762206
# SRR17762207
# SRR17762208

# ######### Monkey ############
# SRR17762209
# SRR17762210
# SRR17762211
# SRR17762212
# SRR17762213
# SRR17762214

# ######### Mouse ############
# SRR17762215
# SRR17762216
# SRR17762217
# SRR17762218
# SRR17762219
# SRR17762220
# SRR17762221
# SRR17762222
# SRR17762223

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## load SRA Toolkir
module load sratoolkit

## List current modules for reproducibility
module list

## Locate file
SAMPLE=$(awk "NR==${SGE_TASK_ID}" SRR_Acc_List.txt)
echo "Processing sample ${SAMPLE}"
date

## fastq-dump
# --split-files to split paired reads
# --include-technical to include barcodes
# -O to output directory
fastq-dump --split-files /dcs04/lieber/marmaypag/Yu_crossAMY_LIBD4125/raw-data/SRA/${SAMPLE} -O /dcs04/lieber/marmaypag/Yu_crossAMY_LIBD4125/raw-data/FASTQ

# example of concatenate multiple SRAs for 1 GEO
#cat SRR1234.fastq.gz  SRR1235.fastq.gz  SRR1236.fastq.gz > GSM5678.fastq.gz

# move fastq files into their own folder
mkdir /dcs04/lieber/marmaypag/Yu_crossAMY_LIBD4125/raw-data/FASTQ/${SAMPLE}
mv ${SAMPLE}_1.fastq dcs04/lieber/marmaypag/Yu_crossAMY_LIBD4125/raw-data/FASTQ/${SAMPLE}
mv ${SAMPLE}_2.fastq dcs04/lieber/marmaypag/Yu_crossAMY_LIBD4125/raw-data/FASTQ/${SAMPLE}

## End
echo "**** Job ends ****"
date