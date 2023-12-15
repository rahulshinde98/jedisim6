import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import seaborn as sns
from scipy.integrate import cumtrapz
from scipy.integrate import trapz
from scipy import interpolate

import sys



#=======================
# Note most clusters are at mass scale of 1e14 solar mass
mass_scale = 1.e14
scale_tag = "14"

#-----------------------
if len(sys.argv)!=3:
    print("Usage: python this.py mass_filename hist_figurename")
    sys.exit(1)

mass_filename = sys.argv[1]
hist_figurename = sys.argv[2]

#-----------------------
tag = mass_filename.split('/')[-1].split('.')[-2]

data = np.loadtxt(mass_filename)
data /= mass_scale

plt.figure()
plt.hist(data, bins=20, histtype="step", density=True)
plt.xlabel("M200 [$10^{%s}M_\odot$]"%scale_tag)
plt.ylabel("density")

ax = sns.kdeplot(data, shade=False)
x, y = ax.get_lines()[0].get_data()

# Initial probability is 0

density = interpolate.interp1d(
                            x, 
                            y, 
                            bounds_error=False,
                            fill_value=0.,
                                ) 

cdf_inverse = interpolate.interp1d(
                            cumtrapz(density(x), x, initial=0),
                            x,
                                )



#-----------------------
x_median = cdf_inverse(0.5)
y_median = density(x_median)

x_negative_1sigma = cdf_inverse(0.5 - 0.341)
x_positive_1sigma = cdf_inverse(0.5 + 0.341)

unity = trapz(density(x), x)
x_mean = trapz(x*density(x), x)
y_mean = density(x_mean)

#nearest_half = np.abs(cdf-0.5).argmin()
#nearest_negative_1sigma = np.abs(cdf-(0.5-0.341)).argmin()
#nearest_positive_1sigma = np.abs(cdf-(0.5+0.341)).argmin()
#
#x_median = x[nearest_half]
#y_median = y[nearest_half]
#
#x_negative_1sigma = x[nearest_negative_1sigma]
#y_negative_1sigma = y[nearest_negative_1sigma]
#
#x_positive_1sigma = x[nearest_positive_1sigma]
#y_positive_1sigma = y[nearest_positive_1sigma]
#
#print("y, x: ", y, x)
#unity = trapz(y, x)
#x_mean = trapz(y*x, x)
#print("x_mean: ", x_mean)
#nearest_x_mean = np.abs(x-x_mean).argmin()
#x_mean = x[nearest_x_mean]
#y_mean = y[nearest_x_mean]
print("unity, x_mean, y_mean: ", unity, x_mean, y_mean)
#
plt.vlines(x_median, 0, y_median, linestyles="dashed", label="median", colors="k")
plt.vlines(x_mean, 0, y_mean, linestyles="dotted", label="mean", colors="k")
##plt.vlines(x_negative_1sigma, 0, y_negative_1sigma, linestyles="dashed", label="$-1\sigma$")
##plt.vlines(x_positive_1sigma, 0, y_positive_1sigma, linestyles="dashed", label="$+1\sigma$")
plt.legend()
#
##plt.minorticks_on()

print("x_median: %.3e"%x_median)
print("x_negative_1sigma: %.3e"%x_negative_1sigma)
print("x_positive_1sigma: %.3e"%x_positive_1sigma)

# Note np.sqrt(2.) is for bootstrap
# Reason: Using the full sample would make the errorbars narrower
plt.title(r"%s: $%.2f_{%.2f}^{+ %.2f} \times 10^{%s}M_\odot$"%(
                    tag, 
                    x_median, 
                    (x_negative_1sigma - x_median) / np.sqrt(2.), 
                    (x_positive_1sigma - x_median) / np.sqrt(2.),
                    scale_tag,
                )
        )

plt.savefig(hist_figurename)

np.savetxt(hist_figurename.replace(".png", ".txt"), np.array([x_median, (x_negative_1sigma - x_median) / np.sqrt(2.), (x_positive_1sigma - x_median) / np.sqrt(2.)] )  )
