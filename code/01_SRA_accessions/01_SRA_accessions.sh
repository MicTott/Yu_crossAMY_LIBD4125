#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=150G
#$ -N Yu_crossAMY_LIBD4125_SRA
#$ -o ./logs/SRA.$TASK_ID.txt
#$ -e ./logs/SRA.$TASK_ID.txt
#$ -m e
#$ -t 1-29
#$ -tc 10


# GEO location
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195445

# ######### Human ############
# GSM5836863	AMH1 snRNA-seq
# GSM5836864	AMH2 snRNA-seq
# GSM5836865	AMH3 snRNA-seq

# ######### Monkey ############
# GSM5836866	AMMA1 snRNA-seq
# GSM5836867	AMMA2 snRNA-seq
# GSM5836868	AMMA3 snRNA-seq

# ######### Mouse ############
# GSM5836869	AMMS1  snRNA-seq
# GSM5836870	AMMS2  snRNA-seq
# GSM5836871	AMMS3  snRNA-seq
# GSM5836872	AMMS4  snRNA-seq
# GSM5836873	AMCEA1 snRNA-seq
# GSM5836874	AMCEA2 snRNA-seq
# GSM5836875	AMBLA1 snRNA-seq
# GSM5836876	AMBLA2 snRNA-seq
# GSM5836877	AMBLA3 snRNA-seq

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

## Prefetch 
prefetch ${SAMPLE} --max-size 150GB -O /dcs04/lieber/marmaypag/Yu_crossAMY_LIBD4125/raw-data/SRA

## End
echo "**** Job ends ****"
date