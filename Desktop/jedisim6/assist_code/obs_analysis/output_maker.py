#Sam Ferraro 6/26/23 10:41:30
#this script turns the fiat file into a csv with the correct columns.

#import libraries
import pandas as pd
import sys

#get the prefix
prefix = sys.argv[1]

#get the file
file_name = prefix + "_LSST_convolved_noise_header_readout_filtered.fiat"
fiat = pd.read_csv(file_name, skiprows=11, sep=' ', names=['x','y','dist','eT','eX','e1','e2','sflux','gflux'])

#add a z column
z = [0.9] * fiat.shape[0]
fiat['z'] = z

#export to csv
output_path = "../../../outputs/" + prefix + "/final_" + prefix + ".csv"
fiat.to_csv(output_path,columns=['x','y','e1','e2','z'],index=False)
