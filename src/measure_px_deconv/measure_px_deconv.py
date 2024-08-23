#Reads in the best fit, creates a new Galfit file, instead of having it free it has it fixed, and it doesn't have the PSF in it. 
#Then it measures the size. It measures the size by summing the flux in the model it divides it by two which is like the half light 
#value of the galaxy and then it takes the model image it sorts the pixels in brightness and then it adds the brightness until it 
#equals a half bright value. It's finding which pixels add up to half of the light.

#Written by Anthony Pahl, 7/1/21

## The majority of these modules are probably not needed
import os, subprocess, sys
import shutil
import glob
from copy import copy

import numpy as np
import pandas as pd
from pandas import Series, DataFrame
from astropy.wcs import WCS
from astropy.io import fits, ascii
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.nddata import Cutout2D
import regions

import numpy.ma as ma
import json

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.backends.backend_pdf import PdfPages
from astropy.visualization import SqrtStretch
from astropy.visualization import ZScaleInterval
from astropy.visualization import make_lupton_rgb

from astropy.cosmology import FlatLambdaCDM

from vars import *

#I did research on this and it said it makes a FlatLambda Cold Darm Matter model.  <-- What does this actually mean?

#Hubble Constant: Represents the expansion of the universe.

#Om0=0.3 This sets the matter density parameter to 0.3 it represents the fraction of the universe's total energy density.
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

# This should point to your galfit installation, just like in src/galfit/prep_galfit.py

#To the file or the directory where the installation is?
def galfit(fin):
    cmd = '/Users/rschmidt-eder/Galfit/galfit3.0.7d/galfit {}'.format(fin)  # Updated to point to the actual executable
    #Subprocess is a module that allows you to create new processes in python.

    #A process includes program code that the CPU executes. It has current activity including registers and variables. Each process has a space in memory. It also includes file descriptors and handles for execution.

    #Popen allows more flexibility toward operating a process.

    #cmd is the thing you're executing (in this case Galfit.)

    #It is operating through the shell, allowing for shell commands and versatility.
    sp = subprocess.Popen(cmd, shell=True)

    #Ensures that the script does not proceed until the process completes.
    sp.wait()

    #Returns the directory to Galfit with a new line.
    return cmd + '\n'

# Reads in pandas dataframe for the object information.
pcord = pd.read_pickle(pdata_path + '/field_mem.pkl')
# This pandas dataframe we don't have for JADES, as everything is one field. can probably comment out

#pfield = pd.read_pickle(pdata_path + '/stamps/field_mem.pkl')
#pcord = pcord.merge(pfield, left_index=True, right_index=True)
#pcord['ID_3DHST'] = np.nan
#pcord['F115W_SEx_ID'] = np.nan
#pcord['F277W_SEx_ID'] = np.nan
#pcord['F356W_SEx_ID'] = np.nan
#pcord['F444W_SEx_ID'] = np.nan

#pfits = pd.read_pickle('/home/apahl/hst_lyc_sfrd/scripts/size_err/bestfit_params.pkl')

# These are the parameters we will be calculating. Effective radius (RE), number of pixels (N), half-light pixel value (HL)

#Makes an empty columns list.
cols = []

#Iterates over our list of bands/filters.
for fi in run_strs:
    #For each filter adds a column for effective radii, sersic index, and HL.

    #I.e. "RE_F115W_single"
    cols.append('RE_{}'.format(fi))
    cols.append('N_{}'.format(fi))
    cols.append('HL_{}'.format(fi))

# Create a pandas dataframe to store the results in, the rows for pandas becomes the pickle file index and the columns become what we've appended to columns

#psize becomes a Pandas object.
psize = DataFrame(index=pcord.index, columns=cols)

#Iterows looks like this (kind of like a dictionary with indexes):
#Index: 0
#Row data:
#A    1
#B    4

# For every object in this new iterrows
for obj,r in psize.iterrows():
    # Select object from info dataframe

    #The current position is equal to the location of the current index in the rowed pickle file.
    pcur = pcord.loc[obj]

    #For every filter in the filter list.
    for run in run_strs:

        # The fit file is equal to the processed data path/galfit/???

        #^ Where are the other format arguments?
        fit_file = '{}/galfit/{}/{}/{}.{}.fit.fits'.format(pdata_path,
                                                           obj, run, obj, run)
        # If it doesn't exist (so, if a fit hasn't been generating), skip
        if not os.path.isfile(fit_file):
            continue

        # Read in science, model, and header info from the .FITS file
        try:  # Added try-except block to handle missing FITS files
            hdu = fits.open(fit_file)
        except FileNotFoundError:  # Skip if the file is not found
            print(f"File not found: {fit_file}, skipping.")
            continue  # Skip this iteration

        headsci = hdu[1].header
        head = hdu[2].header
        model = hdu[2].data
        
        # Figure out number of components that were fit for this object
        n_comp = 1
        d_flag = True
        while (d_flag):
            try:
                if 'COMP_{}'.format(n_comp) in head:
                    n_comp += 1
                else:
                    n_comp -= 1
                    d_flag = False
                    break
            except KeyError:
                print('failed on {}'.format(n_comp))
                n_comp -= 1
                break
        if 'sky' in head['COMP_{}'.format(n_comp)]:
            n_comp -= 1

        # Read in sky value, if it exists
        if 'COMP_{}'.format(n_comp + 1) in head:
            skyval = float(head['{}_SKY'.format(n_comp + 1)].split(' ')[0].replace('*', '').replace('[', '').replace(']', ''))
        else:
            skyval = 0
        
        ##re-generate a GALFIT model based on the best fit parameters, but without a PSF
        #find the *.galfit file that corresponds to this best fit. This logfile will have the best fit params in it
        orig_logfile = '{}/galfit/{}/{}/{}'.format(pdata_path, obj, run, head['LOGFILE'])
        with open(orig_logfile,'r') as f:
            config_lines = f.readlines()

        #Change output filename from *.fits to *.dmodel.fits (named for "deconvolved model")
        old_fit_out = os.path.basename(fit_file)
        new_fit_out = old_fit_out.replace('.fits','.dmodel.fits')
        config_lines = [subs.replace(old_fit_out, new_fit_out) for subs in config_lines]

        #Remove the PSF from the input file
        old_psf = './{}.{}.psf.fits'.format(obj, run)
        config_lines = [subs.replace(old_psf, 'none') for subs in config_lines]
        
        #Make all free paramters fixed rather than free (turn the 1s into 0s)
        ix_1 = 39
        ix_n = 13
        #xy
        ix_o = 1
        for ic in range(n_comp):
            ix_l = ix_1 + ix_n*ic + ix_o
            items = config_lines[ix_l].split()
            items[3] = '0'
            items[4] = '0'
            config_lines[ix_l] =' '+' '.join(items) + '\n'
        #rest
        for ic in range(n_comp):
            for ix_o in np.arange(8)+2:
                ix_l = ix_1 + ix_n*ic + ix_o
                items = config_lines[ix_l].split()
                items[2] = '0'
                config_lines[ix_l] =' '+' '.join(items) + '\n'
            if 'psf' in head['COMP_{}'.format(ic+1)]:
                ix_l = ix_1 + ix_n*ic
                items = config_lines[ix_l].split()
                items[1] = 'sersic'
                config_lines[ix_l] =' '+' '.join(items) + '\n'
                ix_l = ix_1 + ix_n*ic + 3
                items = config_lines[ix_l].split()
                items[1] = str(0.2)
                config_lines[ix_l] =' '+' '.join(items) + '\n'
                
        #sky, free to fixed
        if 'COMP_{}'.format(n_comp+1) in head:
            ic = n_comp
            ix_o = 1
            ix_l = ix_1 + ix_n*ic + ix_o
            items = config_lines[ix_l].split()
            items[2] = '0'
            config_lines[ix_l] =' '+' '.join(items) + '\n'
            ix_l = ix_1 + ix_n*ic + ix_o + 3
            items = config_lines[ix_l].split()
            items[1] = '1'
            config_lines[ix_l] =' '+' '.join(items) + '\n'
        
        #write out this galfit input file (*.dmodel.galfit)
        new_config_file = fit_file.replace('.fits','.dmodel.galfit')
        with open(new_config_file, 'w') as f:
            f.writelines(config_lines)

        #run galfit on this new input file
        init_dir = os.getcwd() + '/'
        id_dir = pdata_path + '/galfit/{}'.format(obj)
        run_dir = id_dir + '/{}'.format(run)
        os.chdir(run_dir)
        arg = '{}.{}.fit.dmodel.galfit'.format(obj, run)
        res = galfit(arg)
        print(res)
        os.chdir(init_dir)

        #read in new model fit (deconvolved)
        try:  # Added try-except block to handle missing new FITS file
            hdu = fits.open('{}/{}'.format(run_dir, new_fit_out))
        except FileNotFoundError:  # Skip if the file is not found
            print(f"File not found: {run_dir}/{new_fit_out}, skipping.")
            continue  # Skip this iteration
            
        headsci = hdu[1].header
        head = hdu[2].header
        model = hdu[2].data

        ##subtract sky value
        model = model - skyval
        # calculate the total amount of flux in the entire image
        total_cts = model.sum()
        # we will count the number of pixels (from brightest to dimmest) down this "half-light" value
        hl_val = total_cts/2.

        # calculate the number of pixels 
        mod_vals = model.flatten()
        mod_vals.sort()
        mod_vals = np.flip(mod_vals)
        #integrate up to hl_val
        msum = 0.0
        for i,e in enumerate(mod_vals):
            msum += e
            if msum > hl_val:
                break
        n_pix = i+1
        fi = run.split('_')[0]
        field = pcord.loc[obj, ['field_{}'.format(fi)]].values[0]
        # Convert from number of pixels to square arcsecs
        pxsc = pxsc_dir[fi][field]
        size_sqarcsec = n_pix * pxsc * pxsc
        # Convert from total "effective area" to "effective radius"
        Re = np.sqrt(size_sqarcsec / np.pi)

        # Place measurements into dataframe
        psize.loc[obj, 'RE_{}'.format(run)] = Re
        psize.loc[obj, 'N_{}'.format(run)] = n_pix
        psize.loc[obj, 'HL_{}'.format(run)] = e + skyval

        # If this is an object with only one sersic profile, compare R_e from this method to the simpler method
        if n_comp == 1:
            print(run)
            r_e = float(head['1_RE'].split(' ')[0].replace('*','').replace('[','').replace(']',''))
            q = float(head['1_AR'.format(i)].split(' ')[0].replace('*','').replace('[','').replace(']',''))
            size = np.pi * r_e * r_e * q
            print('Size from galfit: {}px2'.format(size))
            print('Size from isophote: {}px2'.format(n_pix))
        
# Write out pandas dataframe
psize.to_pickle('{}/size/isophot_size_deconv.pkl'.format(res_path))

drops = [col for col in psize.columns if ('sm' in col)]
psizea = psize.drop(drops, axis=1)