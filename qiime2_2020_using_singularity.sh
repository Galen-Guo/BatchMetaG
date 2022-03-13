#!/bin/bash
#SBATCH -c 5                               # Number of CPUS requested. If omitted, the default is 1 CPU.
#SBATCH --mem=50gb                         # mem in gb
#SBATCH -t 4-0:0:0                         # How long will your job run for? If omitted, the default is 3 hours.
#SBATCH -J qiime2                 	# Name of job

cd $HOME

module load singularity
singularity build qiime2-2020.sif docker://qiime2/core:2020.2