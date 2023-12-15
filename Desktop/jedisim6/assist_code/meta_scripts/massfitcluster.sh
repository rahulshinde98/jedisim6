#!/bin/bash
#********************************
# MASS_FIT_GLOBAL JOB INFORMATION
#********************************

# Walltime requested for job (3 hrs)
#SBATCH -t 3:00:00

# Request use of 1 core on a single node (-N 1)
#SBATCH -N 1
#SBATCH --cpus-per-task=8

# Request 16 GB of memory
#SBATCH --mem-per-cpu=4G

# Define Oscar partition to use
#SBATCH -p batch

# mass_fit_global output and error file names
# Use '%x' for Job Name, '%A' for array-job ID, '%j' for job ID and '%a' for task ID`
#SBATCH -e ./slurm_output/%x-%j.err
#SBATCH -o ./slurm_output/%x-%j.out

# Notify user if job fails of ends
#SBATCH --mail-user=samuel_ferraro@brown.edu
#SBATCH --mail-type=FAIL,END

#*****************
# EXECUTE COMMANDS
#*****************
# get time at beginning of execution step
t0=$(date +%s)
# load desired PYTHON module
module load python/3.7.4
source ~/massfit/bin/activate
# commands to be executed
python3 ./assist_code/mass_fit_global/mass_fit_global.py ./outputs/${1}/final_${1}.csv ./outputs/${1}/${1}_massloc.txt fake ./outputs/${1}/${1}_masses.txt 8

#python3 ./assist_code/mass_fit_global/mass_fit_global.py ./outputs/90_${1}/final_90_${1}.csv ./outputs/90_${1}/90_${1}_massloc.txt fake ./outputs/90_${1}/90_${1}_masses.txt 8
