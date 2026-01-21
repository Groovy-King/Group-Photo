#!/bin/bash -l
#SBATCH --job-name=Mass_Estimator
#SBATCH --account=rrg-babul-ad # adjust this to match the accounting group you are using to submit jobs
#SBATCH --time=0-03:00:00         # adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1      # adjust this if you are using parallel commands
#SBATCH --mem=12000             # adjust this according to the memory requirement per node you need

# Choose a version of MATLAB by loading a module:
module load matlab/2025b.1
matlab -singleCompThread -batch "cd('Rafee_Codes/'); Find_PMlG_virial_All(1, 1)"
