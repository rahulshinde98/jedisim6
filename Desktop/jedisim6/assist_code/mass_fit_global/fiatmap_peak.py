# Locate the fiatmap peak
# Return physical pixel coordinate
#   - Get ra, dec of the peak value position on fiatmap by fiatmap WCS
#   - Convert to physical pixel coordinate by coadd WCS
#   - Note HSM origin is (0,0) but fits image is (1,1) so the above coordinate should minus one
# Usage: python this.py fiatmap.fits coadd.fits 
# Output: x, y pixel value and ra, dec, peak fiatmap value


#=======================
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 

from astropy.wcs import WCS 
import astropy.io.fits as pyfits

import sys 



#=======================
if len(sys.argv) != 4:
    print("Usage: python this.py fiatmap_wcs_filename DATA_path fiatmap_peak_pixel_filename")
    sys.exit(1)

fiatmap_filename = sys.argv[1]
#coadd_filename = sys.argv[2]
DATA_path = sys.argv[2]
#patch = sys.argv[3]
fiatmap_peak_pixel_filename = sys.argv[3]

patch_ref = "5,5"
band_ref = 'r'
coadd_filename = "%s/rerun/coadd2/deepCoadd/%s/0/%s.fits"%(DATA_path, band_ref, patch_ref)


#-----------------------
with pyfits.open(fiatmap_filename) as hdul: 
    wcs_fiatmap = WCS(hdul[0].header)

with pyfits.open(coadd_filename) as hdul: 
    wcs_coadd = WCS(hdul[1].header)


#-----------------------
fiatmap = pyfits.getdata(fiatmap_filename)
#plt.imshow(fiatmap)
#plt.savefig("test.png")
row_peak, col_peak = np.unravel_index(fiatmap.argmax(), fiatmap.shape) 
print("row_peak, col_peak: ", row_peak, col_peak)
x_peak = col_peak
y_peak = row_peak
print("x_peak, y_peak: ", x_peak, y_peak)


ra_peak, dec_peak = wcs_fiatmap.all_pix2world([x_peak], [y_peak], 0)
ra_peak = ra_peak[0]
dec_peak = dec_peak[0]
print("ra_peak, dec_peak: ", ra_peak, dec_peak)

#-----------------------
x_peak_new, y_peak_new = wcs_coadd.all_world2pix([ra_peak], [dec_peak], 0)
x_peak_new = x_peak_new[0]
y_peak_new = y_peak_new[0]
x_peak_new += int(patch_ref.split(',')[0])*4000. - 100.
y_peak_new += int(patch_ref.split(',')[1])*4000. - 100.
print("x_peak_new, y_peak_new: ", x_peak_new, y_peak_new)

#-----------------------
np.savetxt(fiatmap_peak_pixel_filename, 
            np.array([
                        x_peak_new, 
                        y_peak_new
                    ])
            )
