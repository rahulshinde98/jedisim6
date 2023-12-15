#!/bin/bash


#this script should setup lsst and then run the dark matter measurement

current_path=$PWD

module load anaconda/3-5.2.0
source activate ~/anaconda/lsst2
source eups-setups.sh
setup lsst_distrib


# Send this to the obs_file path, you may have to change the path
#cd /oscar/data/idellant/sferrar2/obs_file
cd ../../obs_analysis/obs_file
setup -k -r .
scons

cd ${current_path}

source measure_matter.sh $1

echo "ran measure_matter.sh"
