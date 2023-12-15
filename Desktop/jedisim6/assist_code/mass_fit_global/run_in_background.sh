#!/bin/bash
#********************************
# MASS_FIT_GLOBAL JOB INFORMATION
#********************************
# mass_fit_global Job Name
#SBATCH -J MFG_SF

# Walltime requested for job (3 hrs)
#SBATCH -t 3:00:00

# Request use of 1 core on a single node (-N 1)
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1

# Request 16 GB of memory
#SBATCH --mem-per-cpu=4G

# Define Oscar partition to use
#SBATCH -p batch

# mass_fit_global output and error file names
# Use '%x' for Job Name, '%A' for array-job ID, '%j' for job ID and '%a' for task ID`
#SBATCH -e %x-%j.err
#SBATCH -o %x-%j.out

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
python3 mass_fit_global.py /users/sferrar2/data/sferrar2/lensing_simulation_scripts/outputs/bg180000_z09_${1}/final_bg180000_z09_${1}.csv /users/sferrar2/data/sferrar2/lensing_simulation_scripts/outputs/bg180000_z09_${1}/bg180000_z09_${1}_massloc.txt fake output_${1} 8

#python3 mass_fit_global.py /users/sferrar2/data/sferrar2/lensing_simulation_scripts/outputs/90_bg180000_z09_${1}/final_90_bg180000_z09_${1}.csv /users/sferrar2/data/sferrar2/lensing_simulation_scripts/outputs/90_bg180000_z09_${1}/90_bg180000_z09_${1}_massloc.txt fake output_${1} 8
