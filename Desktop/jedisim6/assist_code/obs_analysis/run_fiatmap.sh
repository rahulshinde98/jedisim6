# Variables
# Remember to adjust cat_name, and maybe r_in
image_name=$1
cat_name="${image_name}LSST_convolved_noise_header_readout_filtered.fiat" 
box_width="10"
r_in="100" #"1000"
r_out="3000"
out_fits="${cat_name}_rin${r_in}_rout${r_out}_all.fits"
fiat_map="/oscar/data/idellant/Clusters/fiatmap"
weight_file="/oscar/data/idellant/Clusters/area.balsa"

# Script
${fiat_map} -b ${box_width} -c "x y e1 e2" -w ${weight_file} ${cat_name} ${r_in} ${r_out} ${out_fits}

