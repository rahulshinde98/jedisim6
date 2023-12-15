#!/bin/bash            
#SBATCH --time=10:00:00
#SBATCH -n 1
#SBATCH --mem=4G
#Python 3 Environment                                                                                                                       
module load python/3.7.4

image_name=$1

#extract data from fits table                                                                         
python ../../obs_analysis/readout_src_cal_flux.py output/src/${image_name}LSST_convolved_noise_header/src.fits > ${image_name}LSST_convolved_noise_header_readout.txt

#filter out erraneous lines, might want to check how this affects things/why there are errors?
module load perl
~/data/Clusters/fiat/bin/fiatfilter "e1<2 && e1>-2 && e2<2 && e2>-2 && x>0 && y>0" ${image_name}LSST_convolved_noise_header_readout.txt > ${image_name}LSST_convolved_noise_header_readout_filtered.fiat

#create dark matter map using fiatmap
cp ../../obs_analysis/run_fiatmap.sh run_fiatmap.sh
bash run_fiatmap.sh ${image_name}

#use source extractor to get the position of the measured dark matter 
cp ../../obs_analysis/default.param default.param
cp ../../obs_analysis/default.sex default.sex
cp ../../obs_analysis/default.conv default.conv

module load sextractor
sex ${image_name}LSST_convolved_noise_header_readout_filtered.fiat_rin100_rout3000_all.fits

cd ..
