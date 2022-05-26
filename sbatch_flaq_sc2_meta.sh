#!/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=flaq_sc2_meta
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ENTER EMAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=100gb
#SBATCH --time=3-00
#SBATCH --output=flaq_sc2_meta.%j.out
#SBATCH --error=flaq_sc2_meta.%j.err

#Run script/command and use $SLURM_CPUS_ON_NODE

primers="4.1"

python flaq_sc2_meta.py fastqs_ww/ --primer_bed primers/ARTIC-V${primers}.bed --lib_frag frag --threads $SLURM_CPUS_ON_NODE --ref_fasta reference/nCoV-2019.reference.fasta --ref_gff reference/GCF_009858895.2_ASM985889v3_genomic.gff 
