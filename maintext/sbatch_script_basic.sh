#!/bin/bash
#SBATCH -J warm_1
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --ntasks-per-core=1
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<gaurav.baruah@uni-bielefeld.de>

Rscript cluster_run_warming.R

