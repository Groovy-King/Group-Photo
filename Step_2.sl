#!/bin/bash -l
#SBATCH --job-name=Mass_Estimator
#SBATCH --account=rrg-babul-ad # adjust this to match the accounting group you are using to submit jobs
#SBATCH --time=0-06:00:00         # adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1      # adjust this if you are using parallel commands
#SBATCH --mem=12000             # adjust this according to the memory requirement per node you need
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Choose a version of MATLAB by loading a module:
module load matlab/2025b.1

vol=$1
slice=$2

matlab -singleCompThread -batch "cd('Rafee_Codes/'); Find_PMlG_virial_All(${vol}, ${slice})"
