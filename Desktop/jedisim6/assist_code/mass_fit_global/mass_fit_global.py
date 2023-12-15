# The algorithm comes from Jacqueline McCleary's work
# Goal: globally fit a shear profile with a single-NFW profile
# Usage: python script.py catalog_csv cluster_center_filename cluster_name mass_filename 
# This catalog has columns: x, y [arcsec], g1, g2, z
# Note some structure of the functions was learned from Nan Li's code
# The idea of "minimization" also comes from CLMM



#=======================
import numpy as np
from scipy.optimize import minimize

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt

from astropy.cosmology import FlatwCDM
from astropy import units as u
from astropy.table import Table

import multiprocessing

import sys

#=======================
# Constants
#-----------------------
# Note in python, variables out of functions' definitions are all treated as global
# Here we use big letters to imitate C-language and emphasize they are global variables

#Z_REF = 5.0   # (source) reference redshift: deflection angle depends on redshift
Z_SOURCE = 0.6   # a fixed redshift for background galaxies
RAD_TO_ARCSEC =  206264.80624709636   # 1/1^{''} #arcsec per rad
COSMO = FlatwCDM(H0=71, Om0=0.2648, Ob0=0.04479, w0=-1.000000)
V_C = 2.998e5   # km/s speed of light
G = 4.3011790220362e-09   # Mpc (Msun)^-1 (km/s)^2
ARCSEC_PER_PIX_DECAM = 0.263

BOOTSTRAP_NUM = 1000 #bootstrap number (times)

LOW_BOUND = 14
HIGH_BOUND = 16

#=======================
# Functions
#-----------------------
# Halo

def d_a(z):
    
    '''
        Angular diameter distance [Mpc]
    '''
    
    return COSMO.angular_diameter_distance(z).value


#-----------------------
def d_a2(z1, z2):
    
    '''
        Angular diameter distance: z1<z2
    '''
    
    return COSMO.angular_diameter_distance_z1z2(z1, z2).value


#-----------------------
def rho_crit(z):
    
    '''
        Critical density at z
        Output: unit: M_sun Mpc^-3
    '''
    
    # original unit is g/cm^3
    return COSMO.critical_density(z).to(u.solMass/u.Mpc/u.Mpc/u.Mpc).value


#-----------------------
def r200_M200(m200, z):
    
    '''
        Input: m200: unit M_sun
        Output: unit: Mpc
    '''
    
    return (3.0/4.0*m200/(200.0*rho_crit(z))/np.pi)**(1.0/3.0)


#-----------------------
def c_M200(m200, z):   
    
    '''
        c = r200/rs for NFW
        Child et al. 2018, Table 2
        
        Input: m200: unit M_sun (not 1e14 M_sun !)
        Note: the values used here (for individual and relaxed halos) need to be double checked in the future! May use individual all halos
    '''
    
    #a, d, m = 68.4, -0.347, -0.083 
    a, d, m = 75.4, -0.422, -0.089
    return a*(1.0+z)**d*(m200**m)


#-----------------------
# Lensing

#def alpha_factor_z_ref_to_z_s(z_l, z_ref, z_s): 
#    
#    '''
#        Rescale alpha from z_ref to z_s, 
#        (reference redshift to true source redshift)
#        because shear, kappa, alpha all have "1/Sigma_crit", 
#        which prop to Dl*Dls/Ds
#        Input: 
#            alpha_ref: reference redshift 
#            z_l: lens redshift
#            z_s: source redshift 
#        Output: 
#            The "factor" that alpha for z_ref should be multiplied by
#            (D_ref/D_l_ref*D_l_s/D_s)
#    '''
#    
#    return d_a(z_ref)/d_a2(z_l, z_ref)*d_a2(z_l, z_s)/d_a(z_s)


#-----------------------
def Sigma_crit(z_l, z_s):
    
    '''
        Critical surface density for the case of lens plane at z_l and source plane at z_s
        Output: unit: Msun Mpc^-2
    '''
    
    return V_C**2/(4.*np.pi*G)*d_a(z_s)/(d_a(z_l)*d_a2(z_l, z_s))


#-----------------------
def NFW_true_g_xy(x, y, c, m200, z_l, z_s, x_c=0., y_c=0.):
    
    '''
        x, y: unit arcsec
        Assume: x_c=0., y_c=0. [arcsec]
        Reference: Wright and Brainerd 2000
        Return the true *reduced* shear gamma1,2/(1-kappa)
    '''
    
    r200 = r200_M200(m200, z_l)
    r_s = r200/c
    
    # Small angle approximation
    r = np.sqrt((x-x_c)**2 + (y-y_c)**2)
    # For NFW, define x = r/rs
    xx = np.abs(r/RAD_TO_ARCSEC*d_a(z_l)/r_s)
    
    delta_c = 200.0/3.0*c**3.0/(np.log(1.0+c)-c/(1.0+c))
    rho_c = rho_crit(z_l)
    Sigma_c = Sigma_crit(z_l, z_s)
    
    # A calculation range
    Sigma = np.zeros_like(xx)
    g_tmp = np.zeros_like(xx)
    
    # idx: index
    # 0<X<1
    idxa = (xx>0.0)
    idxb = (xx<1.0)
    idx1 = idxa&idxb
    x1 = 1.0/(xx[idx1]*xx[idx1]-1.0)
    x2_up = 2.0*np.arctanh(np.sqrt((1.0-xx[idx1])/(1.0+xx[idx1])))
    x2_down = np.sqrt(1.0-xx[idx1]*xx[idx1])
    x2 = x2_up/x2_down
    x3 = np.log(xx[idx1]/2.0)
    x4 = 4.0/(xx[idx1]*xx[idx1])
    Sigma[idx1] = 2.0*r_s*delta_c*rho_c*x1*(1.0-x2)
    g_tmp[idx1] = x2*x4 + x4*x3 - 2.0*x1 + 2.0*x2*x1
    
    # X=1
    idx2 = (xx==1.0)
    Sigma[idx2] = 2.0*r_s*delta_c*rho_c/3.0
    g_tmp[idx2] = 10.0/3.0 + 4.0*np.log(0.5)

    # X>1
    idx3 = (xx>1.0)
    x1 = 1.0/(xx[idx3]*xx[idx3]-1.0)
    x2_up = 2.0*np.arctan(np.sqrt((xx[idx3]-1.0)/(1.0+xx[idx3])))
    x2_down = np.sqrt(xx[idx3]*xx[idx3]-1.0)
    x2 = x2_up/x2_down
    x3 = np.log(xx[idx3]/2.0)
    x4 = 4.0/(xx[idx3]*xx[idx3])
    Sigma[idx3] = 2.0*r_s*delta_c*rho_c*x1*(1.0-x2)
    g_tmp[idx3] = x2*x4 + x4*x3 - 2.0*x1 + 2.0*x2*x1
    
    kappa = Sigma/Sigma_c
    gamma = r_s*delta_c*rho_c*g_tmp/Sigma_c
    g = gamma/(1.0-kappa)
    
    angle = np.arctan2((y-y_c), (x-x_c))
    #gamma1 = -gamma*np.cos(2.*angle)
    #gamma2 = -gamma*np.sin(2.*angle)
    g1 = -g*np.cos(2.*angle)
    g2 = -g*np.sin(2.*angle)
    

    return g1, g2 #g, kappa, gamma, gamma1, gamma2


#-----------------------
def NFW_true_g_r(r, c, m200, z_l, z_s):
    
    '''
        r: unit arcsec
        Reference: Wright and Brainerd 2000
        Return the true *reduced* shear gamma1,2/(1-kappa)
    '''
    
    r200 = r200_M200(m200, z_l)
    r_s = r200/c
    
    xx = np.abs(r/RAD_TO_ARCSEC*d_a(z_l)/r_s)
    
    delta_c = 200.0/3.0*c**3.0/(np.log(1.0+c)-c/(1.0+c))
    rho_c = rho_crit(z_l)
    Sigma_c = Sigma_crit(z_l, z_s)

    Sigma = np.zeros_like(xx)
    g_tmp = np.zeros_like(xx)
    
    # 0<X<1
    idxa = (xx>0.0)
    idxb = (xx<1.0)
    idx1 = idxa&idxb
    x1 = 1.0/(xx[idx1]*xx[idx1]-1.0)
    x2_up = 2.0*np.arctanh(np.sqrt((1.0-xx[idx1])/(1.0+xx[idx1])))
    x2_down = np.sqrt(1.0-xx[idx1]*xx[idx1])
    x2 = x2_up/x2_down
    x3 = np.log(xx[idx1]/2.0)
    x4 = 4.0/(xx[idx1]*xx[idx1])
    Sigma[idx1] = 2.0*r_s*delta_c*rho_c*x1*(1.0-x2)
    g_tmp[idx1] = x2*x4 + x4*x3 - 2.0*x1 + 2.0*x2*x1
    
    # X=1
    idx2 = (xx==1.0)
    Sigma[idx2] = 2.0*r_s*delta_c*rho_c/3.0
    g_tmp[idx2] = 10.0/3.0 + 4.0*np.log(0.5)

    # X>1
    idx3 = (xx>1.0)
    x1 = 1.0/(xx[idx3]*xx[idx3]-1.0)
    x2_up = 2.0*np.arctan(np.sqrt((xx[idx3]-1.0)/(1.0+xx[idx3])))
    x2_down = np.sqrt(xx[idx3]*xx[idx3]-1.0)
    x2 = x2_up/x2_down
    x3 = np.log(xx[idx3]/2.0)
    x4 = 4.0/(xx[idx3]*xx[idx3])
    Sigma[idx3] = 2.0*r_s*delta_c*rho_c*x1*(1.0-x2)
    g_tmp[idx3] = x2*x4 + x4*x3 - 2.0*x1 + 2.0*x2*x1
    
    kappa = Sigma/Sigma_c
    gamma = r_s*delta_c*rho_c*g_tmp/Sigma_c
    g = gamma/(1.0-kappa)

    return g #, kappa, gamma


#-----------------------
#def kappa_to_mass(x_grid, y_grid, kappa, R, z_l, z_s):
#    '''
#        Convert kappa 2D matrix to enclosed projection mass
#        Assume center is the origin
#        Input: x_grid, y_grid: 2D, arcsec
#            R: arcsec (a number)
#        Output: mass: Msun
#        May be useful for multi-NFW/Cluster-galaxy lensing fit
#    '''
#    Sigma_c = Sigma_crit(z_l, z_s)
#    Sigma = kappa*Sigma_c  # Unit: Msun Mpc^-2
#    r = np.sqrt(x_grid**2+y_grid**2)
#    area_unit = (ARCSEC_PER_PIX_ORIGINAL/RAD_TO_ARCSEC*d_a(z_l))**2
#    mass = np.sum(Sigma[r<R])*area_unit
#    
#    return mass
#
#
##-----------------------
#def kappa_to_g(x_grid, y_grid, kappa, R, dr=1.):
#    '''
#        Input: 
#                kappa 2D matrix
#                R: arcsec (a number)
#                dr=1.0 arcsec
#                x_grid, y_grid: 2D, arcsec
#        Output: g (a number at a given R)
#        May be useful for multi-NFW/Cluster-galaxy lensing fit
#    '''
#    
#    r = np.sqrt(x_grid**2+y_grid**2)
#    kappa_mean = np.mean(kappa[r<R])
#    ind = r<(R+dr)
#    ind &= r>(R-dr)
#    kappa_r = np.mean(kappa[ind])
#    #kappa_r = np.median(kappa[ind])
#    gamma = kappa_mean - kappa_r
#    
#    g = gamma/(1.-kappa_r)
#    
#    return g, kappa_r, gamma


#-----------------------
# For mass fitting
def diff_to_minimize(log10_m200, args):

    """
        Find a difference (positive) for minimization in fitting
    """

    m200 = 10**log10_m200
    
    [x_obs, y_obs, g1_obs, g2_obs, z_l, z_s, x_c, y_c] = args
#    [r_obs, gt_obs, z_l, z_s] = args
    

    #-------------------
    # Fit g1, g2 individually

    x_theory, y_theory = x_obs, y_obs
    # Usage: NFW_true_g_xy(x, y, c, m200, z_l, z_s, x_c=0., y_c=0.)
    g1_theory, g2_theory = NFW_true_g_xy(
                                x_theory, 
                                y_theory, 
                                c_M200(m200, z_l), 
                                m200, 
                                z_l, 
                                z_s,
                                x_c,
                                y_c,
                            )
    
    tmp = (g1_obs - g1_theory)**2 + (g2_obs - g2_theory)**2

#    #-------------------
#    # Fit gt
#    # Usage: NFW_true_g_r(r, c, m200, z_l, z_s)
#    r_obs = np.sqrt( (x_obs - x_c)**2 + (y_obs - y_c)**2 ) 
#    r_theory = r_obs
#    
#    angle = np.arctan2( (y_obs - y_c), (x_obs - x_c) )
#    gt_obs = - g1_obs * np.cos(2.0 * angle) - g2_obs * np.sin(2.0 * angle) 
#    gt_theory = NFW_true_g_r(
#                                r_theory, 
#                                c_M200(m200, z_l), 
#                                m200, 
#                                z_l, 
#                                z_s
#                            )
#    
#    # We assume statistical error of the galaxy sample is much larger than measurement error of single galaxies
#    tmp = (gt_obs - gt_theory)**2 #gt_err_obs**2

    return np.sum(tmp) 



#=======================
# Script 
#-----------------------
# Fitting 
# Note the x,y pixel in DM catalog has been "warped" 
# so they are on a flat plane

if len(sys.argv)!=6:
    print("Usage: python this_script.py catalog_csv cluster_center_filename cluster_name mass_filename cpun")
    sys.exit(1)

csv_name = sys.argv[1]
#data = np.genfromtxt(csv_name, 
#                    dtype=None,
#                    delimiter=',',
#                    names=True)

data = Table.read(csv_name, format="ascii.csv", delimiter=",", names=["x","y","e1","e2","z"])
print(data.info)

#csv_name_tag = csv_name.split('.')[-2]
#mass_hist_figurename = csv_name_tag + "_mass.png"

cluster_center_filename = sys.argv[2]
# In pixel
cluster_center = np.loadtxt(cluster_center_filename)

cluster_name = sys.argv[3]

mass_filename = sys.argv[4]

cpun = int(sys.argv[5])


#-----------------------
# Cluster redshift
z_cl = 0.09
#D_a = d_a(z_cl)

#r = data['r']
#gt = data['gt']
#gt_err = data['gt_err']
#gx = data['gx']
#gx_err = data['gx_err']

# Read columns from csv
x_pixel = data["x"]
y_pixel = data["y"]
e1_obs = data["e1"]
e2_obs = data["e2"]
g1_obs = e1_obs /2
g2_obs = e2_obs /2
z_s = data["z"]
 
x_c_pixel = cluster_center[0]
y_c_pixel = cluster_center[1]

# Convert arcsec (sferrar2: got rid of conversion because I'm sending it in as arcsec)
x_obs = x_pixel * ARCSEC_PER_PIX_DECAM
y_obs = y_pixel * ARCSEC_PER_PIX_DECAM

x_c = x_c_pixel *1400/6144 * ARCSEC_PER_PIX_DECAM
y_c = y_c_pixel *1400/6144 * ARCSEC_PER_PIX_DECAM

print("-------------------------------------------------")
print("Mean of array 'x_obs' elements: ", np.mean(x_obs))
print("Standard Deviation of array 'x_obs' elements: ", np.std(x_obs))
print("Range of array 'x_obs' elements: (" + str(np.min(x_obs)) + ", " + str(np.max(x_obs)) + ")")
print("")
print("Mean of array 'y_obs' elements: ", np.mean(y_obs))
print("Standard Deviation of array 'y_obs' elements: ", np.std(y_obs))
print("Range of array 'y_obs' elements: (" + str(np.min(y_obs)) + ", " + str(np.max(y_obs)) + ")")
print("")
print("-------------------------------------------------")
print("CENTER: (" + str(x_c) + ", " + str(y_c) + ")")
print("Number: ", len(x_obs))

#-----------------------
#r_obs = np.sqrt( (x_obs - x_c)**2 + (y_obs - y_c)**2 )
#angle = np.arctan2( (y_obs - y_c), (x_obs - x_c) )
#gt_obs = - g1_obs * np.cos(2.0 * angle) - g2_obs * np.sin(2.0 * angle)

#-----------------------
print('Fitting...')

test_log10_m200 = 14.
#mass = []
#for test_log10_m200 in range(14, 16):

def run_process(arg):
#    fitting_result = minimize(
#                                diff_to_minimize, 
#                                test_log10_m200, 
#                                args=[r, gt, gt_err, z_cl, Z_SOURCE]
#                            )
    # Bootstrap
    #print("Test 1")
    np.random.seed()
    #print("Test 2")
### Make a list of Falses the length of the data set 
    select = np.zeros(len(x_obs), dtype=bool)
#    select[:50000] = True
    #print("Test 3")
### Set half of them at random to be true
    select_index = np.random.choice(range(len(select)), size=int(len(select)/2) )
    #print("Test 4")
    for index in select_index: select[index] = True
### Grab only those data. Now we have a random sample of half of the data
    x_obs_tmp = x_obs[select]
    y_obs_tmp = y_obs[select]
    g1_obs_tmp = g1_obs[select]
    g2_obs_tmp = g2_obs[select]
    z_s_tmp = z_s[select]
    #print("Test 5")
#    r_obs_tmp = r_obs[select]
#    gt_obs_tmp = gt_obs[select]
    
    #print(type(x_obs_tmp),type(y_obs_tmp),type(x_c),type(y_c))
### Now we minimize the function diff_to_minimize with a starting value test_log10_m200 and bounds
    fitting_result = minimize(
                                diff_to_minimize, 
                                test_log10_m200, 
#                                args=[r_obs_tmp, gt_obs_tmp, z_cl, z_s_tmp],
                                args=[x_obs_tmp, y_obs_tmp, g1_obs_tmp, g2_obs_tmp, z_cl, z_s_tmp, x_c, y_c],
                                bounds=( (LOW_BOUND, HIGH_BOUND), ),
                            )
    

#    print('At test_log10_m200 = %.3e, \
#            log10_fitted mass = %.3e, \
#            fitted mass = %.3e [M_sun]'%(
#                    test_log10_m200, 
#                    fitting_result.x[0],
#                    10**fitting_result.x[0]
#                )
#        )

    # Return mass with unit of solar-mass
    return 10**fitting_result.x[0]

#    mass.append(10**fitting_result.x[0])

#print("mass list: ", mass)
#cpun = 20
p = multiprocessing.Pool(cpun)
result = p.map(run_process, range(BOOTSTRAP_NUM) )
p.close()
p.join()
#print("result: ", result)

'''
# Grab the number
import os
script_path = os.path.abspath(__file__)
script_dir = os.path.split(script_path)[0]
rel_path = "output_files/number.txt"
num_file_path = os.path.join(script_dir, rel_path)
f1 = open(num_file_path, 'r')
num = f1.readlines()[0]
f1.close()
'''
# Write to file
f2 = open(mass_filename, 'w')
for i in range(len(result)):
    f2.write("%e\n"%result[i])
f2.close()
'''
# Increment number
f3 = open(num_file_path, 'w')
f3.write(str(int(num.split()[0]) + 1))
f3.close()
'''
print("DONE")
#m200_theory = 10**fitting_result.x[0]
#gt_theory = NFW_true_g_r(
#                            r, 
#                            c_M200(m200_theory, z_cl), 
#                            m200_theory, 
#                            z_cl, 
#                            Z_SOURCE
#                        )[0]
#
#print("r [arcsec]: ", r)
#print("gt: ", gt)
#print("gt_theory: ", gt_theory)
#
#plt.errorbar(r, gt, gt_err, fmt='-', label='<eT>/2')
#plt.errorbar(r, gx, gx_err, fmt='-', label='<eX>/2')
#plt.plot(r, gt_theory, '--', label='fit:%.1e$M_{\odot}$ zs=%.1f'%(m200_theory, Z_SOURCE))
#
#plt.legend()
#plt.xlabel('CL radius [arcsec]')
#plt.ylabel('<ei>/2')
#plt.minorticks_on()
#
#plt.hlines(0, np.min(r), np.max(r), linestyles=':', colors='gray')
#
#plt.title('%s z=%.3f'%(cluster_name, z_cl))
#
#plt.tight_layout()
#
#plt.savefig(csv_name_tag+'_mass_fit_%s.png'%cluster_name)



#=======================
#-----------------------
# Test
#xdata = np.linspace(100, 1000, 10)
#ydata = NFW_true_g_r(xdata, c=4, m200=1e15, z_l=0.05, z_s=1.0)[0]
#print(xdata)
#print(ydata)
#
#def diff(log10_m200, r):
#    m200 = 10**log10_m200
#    model = NFW_true_g_r(r, 4, m200, z_l=0.05, z_s=1.0)[0]
#    tmp = model-ydata
#    return np.sum(tmp**2) #tmp
#
#for test_log10_m in range(13, 18):
#    a = minimize(diff, test_log10_m, args=xdata)
#    print(a.x[0])
