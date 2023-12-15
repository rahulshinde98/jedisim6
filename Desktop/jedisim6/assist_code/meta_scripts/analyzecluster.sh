#!/bin/bash

# Request an hour of runtime:
#SBATCH --time=2-00:00:00

# Request use of 1 core on a single node (-N 1)
#SBATCH -N 1
#SBATCH -p bigmem

# Specify an output file
#SBATCH -o ./slurm_output/%x-%j.out
#SBATCH -e ./slurm_output/%x-%j.err

# Notify user if job fails of ends
#SBATCH --mail-user=samuel_ferraro@brown.edu
#SBATCH --mail-type=FAIL,END

#============================================#

#Move to the outpur file
cd ./assist_code/jedisim/${1}

#Copy needed files into the output folder
cp ../../obs_analysis/setup_lsst.sh setup_lsst.sh
cp ../../obs_analysis/measure_matter.sh measure_matter.sh
cp ../../obs_analysis/output_maker.py output_maker.py
cp ../../obs_analysis/create_catalog.sh create_catalog.sh

#Analyze
bash setup_lsst.sh ${1}_

#make a csv file/output folder
module load anaconda
bash create_catalog.sh ${1}_

mkdir ../../../outputs/${1}
python output_maker.py ${1}
  
echo 'output made'
cd ..
    
#repeat for the 90_ folder
#cd 90_${1}
#cp ../../obs_analysis/setup_lsst.sh setup_lsst.sh
#cp ../../obs_analysis/measure_matter.sh measure_matter.sh
#cp ../../obs_analysis/output_maker.py output_maker.py
#cp ../../obs_analysis/create_catalog.sh create_catalog.sh
   
#bash setup_lsst.sh 90_${1}_
#module load anaconda
#bash create_catalog.sh 90_${1}_
  
#mkdir ../../outputs/90_${1}
#python output_maker.py 90_${1}

#echo '90_ output made'
#cd ..

#copy important info to output folder
cp ${1}/${1}_catalog.txt ../../outputs/${1}/${1}_catalog.txt
cp ${1}/${1}_LSST_convolved_noise_header.fits ../../outputs/${1}/${1}_noise_header.fits
cp physics_settings/config_z09_${2} ../../outputs/${1}/${1}_config.txt
cp physics_settings/lenses_${2}.txt ../../outputs/${1}/${1}_massloc.txt
  
#cp 90_${1}/90_${1}_catalog.txt ../../outputs/90_${1}/90_${1}_catalog.txt
#cp 90_${1}/90_${1}_convolved_noise_header.fits ../../outputs/90_${1}/90_${1}_noise_header.fits
#cp physics_settings/config_z09_${2}.txt ../../outputs/90_${1}/90_${1}_config.txt
#cp physics_settings/lenses_${2}.txt ../../outputs/90_${1}/90_${1}_massloc.txt

#Clean up everything else
echo "Done, deleting..."
rm -r ${1}
#rm -r 90_${1}
