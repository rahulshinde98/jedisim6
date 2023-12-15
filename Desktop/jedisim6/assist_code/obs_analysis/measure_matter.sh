#!/bin/bash

#this script contains the steps to measure the dark matter
# Assume now you are in the folder that contains simulated images
#must be in lsst environment with obsfile set up

#to run: bash measure_matter.sh image.fits
rm -rf input
mkdir input

echo "lsst.obs.file.FileMapper" > input/_mapper

image_name=$1

ingestImages.py input/ ${image_name}LSST_convolved_noise_header.fits --mode link

> img.list
echo "--id filename=${image_name}LSST_convolved_noise_header.fits" >> img.list

cpun=1
processCcd.py input/ @img.list -j $cpun --timeout 9999999 --configfile ${HOME}/data/config_new --output output
