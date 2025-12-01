#! /bin/bash 

#BSUB -o nextflow.o.%J 
#BSUB -e nextflow.e.%J 
#BSUB -G team301 
#BSUB -n 1
#BSUB -q normal 
#BSUB -R select[mem rusage[mem=10000] span[hosts=1] 
#BSUB -M10000

module load singularityce-4.1.0/python-3.11.6 

singularity remote build slseek.sif slseek.def
