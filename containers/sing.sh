#! /bin/bash

#BSUB -o sing.o.%J
#BSUB -e sing.e.%J
#BSUB -G team301
#BSUB -q small
#BSUB -n 5
#BSUB -R "select[mem>10000] rusage[mem=10000] span[hosts=1]"
#BSUB -M10000

module load singularityce-4.1.0/python-3.11.6 

singularity build --remote slseek.sif slseek.def
