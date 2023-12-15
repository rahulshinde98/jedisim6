#modified @ Sat Dec  2 11:54:14 EST 2017 add child filter
#modified add star filter (psf candidate)

import re
import sys
import numpy as np
from astropy.io import fits as pyfits

if len(sys.argv) != 2:
    print(sys.argv)
    print("len=%d"% len(sys.argv))
    print(sys.stderr, "Usage: python readout_src.py {source_fits_file} > {catalog_file}")
    sys.exit(1);
srcfits = sys.argv[1]

# Load sources and print all columns
if srcfits=='':
    srcfits = 'src.fits'

data_table, header_table = pyfits.getdata('%s' %(srcfits), 1, header=True)

cols = (
        "ext_shapeHSM_HsmSourceMoments_x",
        "ext_shapeHSM_HsmSourceMoments_y",
        "ext_shapeHSM_HsmShapeRegauss_e1",
        "ext_shapeHSM_HsmShapeRegauss_e2",
        "deblend_nChild",    #should be zero
        "base_SdssShape_flux",
        "base_GaussianFlux_flux",
        "flags"
        )

print("# fiat 1.0")
print("# TTYPE1 = x")
print("# TTYPE2 = y")
print("# TTYPE3 = distance")
print("# TTYPE4 = eT")
print("# TTYPE5 = eX")
print("# TTYPE6 = e1")
print("# TTYPE7 = e2")
print("# TTYPE8 = sdss_flux")
print("# TTYPE9 = gauss_flux")

vecs = []
for col in cols:
    v = data_table.field(col)
    vecs.append(v)

zipdata = zip(*vecs)
for vals in zipdata:
    #print(type(vals[4]), vals[4])
    #print(' '.join(map(str, vals)))
    if vals[-1][1]==False and vals[4]==0: 
#not psf candidate, no child
        x, y, e1, e2 = vals[0], vals[1], vals[2], vals[3]
        sdss_flux, gauss_flux = vals[5], vals[6]
        dx, dy = x-2803/2.+0.5, y-2803/2.+0.5
        d = np.sqrt(dx*dx+dy*dy)
        angle = np.arctan2(dy, dx) 
        et = -e1*np.cos(2*angle)-e2*np.sin(2*angle) 
        ex = e1*np.sin(2*angle)-e2*np.cos(2*angle)
        print(' '.join([str(x), str(y), str(d), str(et), str(ex), str(e1), str(e2), str(sdss_flux), str(gauss_flux)]))
