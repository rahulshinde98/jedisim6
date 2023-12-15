#!/bin/bash

# Request an hour of runtime:
#SBATCH --time=2-00:00:00

# Request use of 1 core on a single node (-N 1)
#SBATCH -N 1
#SBATCH --cpus-per-task=15

# Specify an output file
#SBATCH -o ./slurm_output/%x-%j.out
#SBATCH -e ./slurm_output/%x-%j.err

# Notify user if job fails of ends
#SBATCH --mail-user=samuel_ferraro@brown.edu
#SBATCH --mail-type=FAIL,END

#ADDING IN COMMENTS ON THESE TWO
module load cfitsio
module load fftw

#============================================#
cd ./assist_code/jedisim

#Set Variables
galaxies=${8}
npix=$((12288*${1}))
mass=${5}
num=${2}
threads=${9}
noise=${4}

#Generate names to replace in config and lens files.
old_output_folder="bg138000_z09_tmp"
new_output_folder="${3}"
old_prefix="bg138000_z09_tmp_"
new_prefix="${3}_"
old_config="config_z09_tmp"
new_config="config_z09_${num}"
old_gal="num_galaxies=138000"
new_gal="num_galaxies=${galaxies}"
old_lens="lenses_20.000000_4.000000.txt"
new_lens="lenses_${num}.txt"
old_pix="=12288"
new_pix="=${npix}"
old_noise="noise_mean=1"
new_noise="noise_mean=${noise}"
lens="${6} ${7} 2 ${mass} 4"

echo "start bg ${num}......"

#Make config and lens settings for the new run
sed "s/${old_output_folder}/${new_output_folder}/g;s/${old_prefix}/${new_prefix}/g;s/${old_gal}/${new_gal}/g;s/${old_lens}/${new_lens}/g;s/${old_pix}/${new_pix}/g;s/${old_noise}/${old_noise}/g" ./physics_settings/${old_config} > ./physics_settings/${new_config}
echo ${lens} > ./physics_settings/${new_lens}

#Run Jedimaster
echo "jedimaster......"
python -u ./jedimaster.py ./physics_settings/${new_config} ${threads}

echo "bg ${i} finished!"

#python 3 environment
module load python/3.7.4
source ~/massfit/bin/activate

#run addstar.py, add header
cd ./${new_output_folder}
cp ../psf_scalednew_263.fits psf_scalednew_263.fits
python ../addstar.py ${new_prefix}LSST_convolved_noise.fits psf_scalednew_263.fits ${noise}
python ../addwcs2.py ${new_prefix}LSST_convolved_noise.fits ${new_prefix}LSST_convolved_noise_header.fits

cd ..

echo 'Stars added'

#repeat for the 90_ folder
#cd 90_${new_output_folder}
#cp ../psf_scalednew_263.fits psf_scalednew_263.fits
#python ../addstar.py 90_${new_prefix}LSST_convolved_noise.fits psf_scalednew_263.fits
#python ../addwcs2.py 90_${new_prefix}LSST_convolved_noise.fits $90_{new_prefix}LSST_convolved_noise_header.fits

#cd ..

#echo 'Stars added to 90'
    
deactivate
