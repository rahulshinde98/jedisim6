# Lensing Simulation and Analysis README

ssh <username>@ssh.ccv.brown.edu

################################################################################################
# MAKECLUSER
This is how you can run all three parts:
1) Jedisim
2) Obs_file Analysis
3) Mass_fit_global

You will need to do a little set-up for each (see below)
But once they are set up, you can run everything from here.
Call it by running:
     sbatch makecluster.sh <options>

Add whatever options you want to make your image.
To see options call the help option:
   bash makecluster.sh -h
(Bashing it will make it appear in the command line)

The output will be placed in a folder within the outputs directory.
The name of that folder will be bg<number of galaxies>_z09_<number>,
where the final number exists to differentiate files.

################################################################################################
1) JEDISIM

Jedisim generates a fits file simulating the effect of a lensing object
in the foreground on a number of background galaxies. Jedimaster.py is
the control file for this. It calls executables, written in C. (see information
below about these executables)

Jedisim has a number of parameters which can be controlled by the makecluster
file options. You can also access them in ./assist_code/physics_settings/config_z09_tmp.
This file controls aspects from noise, to image size, to the number of
background galaxies. Additionally, generatecluster.sh creates a file in the same
directory called lenses_<number>.txt which determines the location of the lenses
(x,y), the type of lens (currently set to 2 for standard profile), the mass, and
the concentration (currently set to 4).

#===========================#              |SETUP|              #============================#
The old readme claims you need the following packages:
    pip3 install --upgrade astropy astroquery ipython matplotlib numpy pandas requests scipy

If you want, you can make a virtual environment before running that command:
   module load python/3.7.4
   cd ~
   virtualenv run_python3
   source ~/run_python3/bin/activate

I don't use that virtual environment, which implies, I think, I have those packages installed
in my main environment. If you want to use the virtual one, you will have to activate the
virtual environment before jedimaster.py runs in generatecluster.sh. (and deactivate it after)
#============================================================================================#

################################################################################################
1.5) ADDSTAR & ADDWCS
This step adds psf stars to the jedisim image as well as a header.
This is needed for the analysis.

#===========================#              |SETUP|              #============================#
When I was setting this up, I had already set up mass_fit_global so I had that virtual
environment set up already. If you plan to set that up, you don't have to do anything here.

If you aren't going to use massfit, you still need that virtual environment (see massfit set
up for details). If you call it something different, don't forget to change the name in
generatecluster.sh.
#============================================================================================#

################################################################################################
2) OBS_FILE ANALYSIS

This section analyzes the output file of jedisim as though it were an image
being processed by the pipeline. It runs off of an old version of the LSST
pipeline.

In the end, it outputs a file that contains all of the locations and ellipticities
of detected background galaxies.

#===========================#              |SETUP|              #============================#
The set-up for this took me a week. Hopefully, with this guidence, it won't take you as long.

In your home directory, run the following code:
   module load anaconda/3-5.2.0
   conda config --add channels https://conda.lsst.codes/stack/0.13.0
   conda create -name lsst2 python=2
   source activate lsst2
   conda install lsst-distrib

If you run into an error with the time package, run:
   unset EUPS_PATH EUPS_DIR

If you run into an error with astropy and numpy not getting along, it is because even in the
environment, conda is grabbing the new version of astropy from your installed python libraries
You have to change the python file path somewhere in the environment's code. (For this one,
I suggest asking for CCV support and telling them this is the way to fix it. They should be
able to help.)

For your sake, I hope nothing has changed since I did this. If it has, I wish you luck.
(Unless they updated it to the modern version, in which case it should be easier)
#============================================================================================#

################################################################################################
3) MASS_FIT_GLOBAL

This section analyzes the output of the obs_file analysis and attempts
to calculate the mass of the lensing object. This is code from the main
pipeline, so if that has been updated, you might have to update this too.

#===========================#              |SETUP|              #============================#
For this, you just need to set-up a virtual environment. If you don't want to, you can just
install the packages in your main environment. But, you'll have to make an edit in
massfitcluster.sh to remove the call to the virtual environment.

In your home directory run the following code:
   module load python/3.7.4
   cd ~
   virtualenv massfit
   source ~/massfit/bin/activate

Then install the following packages to that virtual environment:
     numpy scipy astropy multiprocessing sys   
#============================================================================================#

################################################################################################
# JEDISIM SOURCE CODE

The jedisim source code is written in C. After my recent changes, there are 6 files:
    i)   jedicatalog_s
    ii)  jedimakeimage_s (jeditransform, jedidistort, jedipaste)
    iii) jediconvlove
    iv)  jedipaste
    v)   jedirescale
    vi)  jedinoise
Which are all called by jedimaster.

i)   The catalog file generates a catalog of background galaxies with random sizes,
     magnitudes, ellipticities and redshifts. To each, it assigns a postage stamp from
     the database. It outputs all of this into a file.
ii)  Makeimage takes that catalog and does three things. It takes the assigned postage
     stamp and transforms it to the set specifications. It then distorts is based on how
     the light would be lensed by the lensing mass. It then pastes it into an array.
     When done, it outputs this array as a fits file. This process used to be 3 parts. In
     between each step, it would output each image as a distinct fits file. While this
     saved runtime memory, it took much longer. I combined them and also added multithreading
     to increase efficiency (at the cost of memory).
iii) Convolve changes the image to make it seem like the telescope took the image.
iv)  Paste puts together the output images from convolve (which takes apart the previous
     image for memory reasons).
v)   Rescale turns the image from HST resolution to LSST resolution. It also chops off
     the borders.
vi)  Noise adds the desired level of noise to the image.

You can read more about the process in jedisim.pdf, which Professor Ian Dell'Antonio can
find for you.

Some of the source code exists in the ./assist_code/jedisim_source/ directory. It can be
compiled by running:
	bash compile.sh <process>
	(ex. bash compile.sh makeimage)

This folder only has the sources I have worked with. The rest can be found at:
     https://github.com/rbliu/jedisim

#===========================#              |SETUP|              #============================#
You only need to do this set-up if you are compiling new jedisim exacutables.

The main thing you need to set up is the cfitsio library. Run this code in the directory where
you are compiling your code:
    git clone https://github.com/healpy/cfitsio.git
    cd cfitsio/
    ./configure --prefix=/users/sferrar2/data/sferrar2/ --disable-curl --enable-reentrant
    make
    make install

That will set up the cfitsio library in a way that should work on CCV and will allow
multithreading. The compile.sh file compiles using the static library of cfitsio. So anyone
should be able to use it, even if they don't install the library themself.
#============================================================================================#


################################################################################################
