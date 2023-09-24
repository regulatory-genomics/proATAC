#!/bin/bash
#SBATCH -J test
#SBATCH -p intel-sc3
#SBATCH -c 2
#SBATCH -q normal
#SBATCH --mem=20G
conda activate  cutadaptenv
/storage/zhangkaiLab/fanjiaqi/data/PUBATAC_pipline/PUMATAC_dependencies/nextflow/nextflow-21.04.3-all -C atac_preprocess.config run proATAC.nf -entry atac_preprocess -resume
