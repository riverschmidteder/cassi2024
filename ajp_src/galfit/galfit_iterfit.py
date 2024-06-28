
import os, subprocess, sys, shutil
import glob
import pandas as pd
from pandas import Series, DataFrame
from astropy.table import Table
import numpy as np
import copy
import multiprocessing as mp
from subprocess import DEVNULL, STDOUT
from getkey import getkey, keys

from astropy.io import fits,ascii
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.nddata import Cutout2D
import astropy.units as u


import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.backends.backend_pdf import PdfPages
from astropy.visualization import SqrtStretch
from astropy.visualization import ZScaleInterval
from astropy.visualization import make_lupton_rgb

from vars import *

sys.path.append("./galfit_iterfit")
from plot_obj import plot_obj
from refit_galfit import refit_galfit
from init_mult_fit import init_mult_fit, remove_mult_fit
from apply_big import apply_big

def find_recentfit(run_dir):
    all_runs = sorted(glob.glob(run_dir + '/galfit.??'))
    return all_runs[-1]

def open_ds9(obj, run):
    fit_file = '{}/galfit/{}/{}/{}.{}.fit.fits'.format(
        pdata_path,
        obj,
        run,
        obj,
        run)

    args = ['ds9', '-frame', 'lock', 'physical','-lock', 'colorbar', 'yes',
            '-cmap', 'sls', 
            fit_file+'[1]', fit_file+'[2]', fit_file+'[3]',
            '-frame', 'first', '-match', 'scalelimits',
            '-zoom', 'to', 'fit']
    #print(' '.join(args))
    subprocess.Popen(args, stdout=DEVNULL)

# Retrieve CEERS objects to fit
pcrs = pd.read_pickle(res_path + '/cat/ceers_vdw.pkl')
# Include both VDW and SExtractor catalog info
pphot = pd.read_pickle(res_path + '/cat/ceers_sex_phot.pkl')
pcrs = pphot.merge(pcrs.iloc[:,4:], on='ID', how='left')
# And self-identified labels
plab = pd.read_pickle(pdata_path + '/galfit/fit_labels.pkl')
plab = plab.dropna(axis=0, how='all')
for row in plab.itertuples():
    obj = int(row.Index.split('.')[0])
    run = row.Index.split('.')[1]
    fi = run.split('_')[0]
    pcrs.loc[obj, 'label_{}'.format(run)] = row.label

# Look for saved file with iterfitting status / comments
if os.path.isfile(pdata_path + '/galfit/iterfit_status.pkl'):
    pst = pd.read_pickle(pdata_path + '/galfit/iterfit_status.pkl')
else:
    pst = Series(False, index=pcrs.index)
if os.path.isfile(pdata_path + '/galfit/iterfit_comment.pkl'):
    pco = pd.read_pickle(pdata_path + '/galfit/iterfit_comment.pkl')
else:
    pco = Series('', index=pcrs.index)
if os.path.isfile(pdata_path + '/galfit/iterfit_complab.pkl'):
    pcl = pd.read_pickle(pdata_path + '/galfit/iterfit_complab.pkl')
else:
    pcl = Series(-1, index=pcrs.index, dtype=object)


# Start by fitting tweak singles, mults, sort by z
# Actually, let's just do everything
###ix_simp = []
###for i,r in pcrs.iterrows():
###    jl = r.label_F115W_single
###    hl = r.label_F160W_single
###    if (jl == 't') or (hl == 't'):
###        ix_simp.append(i)
###ix_mult = []
###for i,r in pcrs.iterrows():
###    jl = r.label_F115W_single
###    hl = r.label_F160W_single
###    if (jl == 'm' or jl == 'r') or (hl == 'm' or hl == 'r'):
###        ix_mult.append(i)
###pcrs = pcrs.loc[set(ix_simp + ix_mult)]
pcrs = pcrs.sort_values(by='z', ascending=False)

###print('To fit with simple tweaks: {}'.format(len(ix_simp)))# - pst[ix_simp].sum()))
###print('To fit with multiz: {}'.format(len(ix_mult)))# - pst[ix_mult].sum())))
print('Total to fit: {}'.format(pcrs.shape[0]))# - pst[pcrs.index].sum()))
ntotal = pcrs.shape[0]# - pst[pcrs.index].sum()


# Open a new emacs window
subprocess.Popen(['emacs'])
input("Start server in emacs (M-x server-start), then press Enter to continue...")
      
    
    


# Define figure
num_x = 5
num_y = len(run_strs)
fig = plt.figure(figsize = (4*num_x+0.5, 4*num_y+1))

for row in pcrs.itertuples():

    # See if already fit
    if pst[row.Index]:
        continue

    # Check how many runs currently exist
    runs = copy.copy(run_strs)
    todel = []
    for i,run in enumerate(run_strs):
        fi = fi_dir[run]
        fit_file = '{}/galfit/{}/{}/{}.{}.fit.fits'.format(
            pdata_path, row.Index, run, row.Index, run)
        if not os.path.isfile(fit_file):
            todel.append(i)
    for i in todel[::-1]:
        del runs[i]
            
    # Plot all fit images
    plot_obj(row.Index, runs, row, fig=fig)
    plt.draw()
    plt.pause(0.001)
    #plt.close(fig)

    # Vars for fitting
    ix_run = 0
    new_run_fl = True
    next_run_fl = False
    quit_fl = False

    while not next_run_fl:
        run = runs[ix_run]
        fi = fi_dir[run]
        field = getattr(row, 'field_{}'.format(fi))
        id_dir = pdata_path + '/galfit/{}'.format(row.Index)
        run_dir = id_dir + '/{}'.format(run)
        init_file = run_dir + '/{}.{}.init.galfit'.format(row.Index, run)
        run_file = run_dir + '/{}.{}.run.galfit'.format(row.Index, run)
        fit_file = run_dir + '/{}.{}.fit.fits'.format(row.Index, run)

        # Skip if label is n
        if 'single' in run:
            if (getattr(row, 'label_{}'.format(run)) == 'n'):# or (getattr(row, 'label_{}'.format(run)) == 's'):
                ix_run+=1
                if ix_run > (len(runs) - 1):
                    next_run_fl = True
                    pst[row.Index] = True
                continue
        elif 'vdw' in run:
            ix_run+=1
            if ix_run > (len(runs) - 1):
                next_run_fl = True
                pst[row.Index] = True
            continue

        if new_run_fl:
            # Open run file in emacs (.run.)
            subprocess.Popen(['emacsclient', run_file], stdout=DEVNULL, stderr=STDOUT)

            # Display best-fit information (include label)
            with fits.open(fit_file) as hdu:
                model_header = hdu[2].header
            mag = model_header['1_MAG'].split(' ')[0]
            if '1_RE' in model_header:
                r_e = model_header['1_RE'].split(' ')[0]
            else:
                r_e = np.nan
            if '1_N' in model_header:
                n = model_header['1_N'].split(' ')[0]
            else:
                n = np.nan
            if '1_PA' in model_header:
                pa = model_header['1_PA'].split(' ')[0]
            else:
                pa = np.nan
            if '1_AR' in model_header:
                q = model_header['1_AR'].split(' ')[0]
            else:
                q = np.nan
            print('===================')
            print(row.Index, '\t {}/{}'.format(pcrs.index.get_loc(row.Index), ntotal))
            print(run+'\n')
            print('Fit Parameters:')
            print('mag={}\nr_e={}\nn={}\npa={}\nq={}\n'.format(mag, r_e, n, pa, q))

        # Ask for input
        inp = input('Action? ')

        if inp == 'f':
            # Refit
            refit_galfit(row.Index, run)
            plt.clf()
            plot_obj(row.Index, runs, row, fig=fig)
            plt.draw()
            plt.pause(0.001)
            new_run_fl = True
        elif inp == 'a':
            # Copy best-fit params to current runfile
            recent_file = find_recentfit(run_dir)
            shutil.copyfile(recent_file, run_file)
            new_run_fl = True
        elif inp == 'r':
            # Copy initial guess to current runfile
            shutil.copyfile(init_file, run_file)
            new_run_fl = True
        elif inp == 'm':
            # Initialize multiple fit
            newrun = init_mult_fit(row.Index, fi, row)
            if not pd.isnull(newrun):
                refit_galfit(row.Index, newrun)
                runs.insert(ix_run+1, newrun)
                plt.clf()
                plot_obj(row.Index, runs, row, fig=fig)
                plt.draw()
                plt.pause(0.001)
            else:
                print('Multiple runs was already initialized')
            new_run_fl = False
        elif inp == 'rm':
            # Remove multiple fit
            delrun = remove_mult_fit(row.Index, fi, row)
            if not pd.isnull(delrun):
                del runs[runs.index(delrun)]
                plt.clf()
                plot_obj(row.Index, runs, row, fig=fig)
                plt.draw()
                plt.pause(0.001)
                if run == delrun:
                    ix_run-=1
            else:
                print('No multiple run found')
            new_run_fl = False
        elif inp == 'bo':
            apply_big(row.Index, fi, row, run)
            new_run_fl = True
        elif inp == 'b':
            # Go back
            ix_run-=1
            if ix_run < 0:
                ix_run = 0
                print('Cant go back any more')
            new_run_fl = True
        elif inp == 'c':
            # Add comment
            comment = input('Enter comment: ')
            pco[row.Index] = comment
            new_run_fl = False
        elif inp == 'co':
            # Delineate which components are part of the target
            complab = input('Enter which components are associated with the target: ')
            pcl[row.Index] = complab
        elif inp == 'd':
            # Open DS9
            open_ds9(row.Index, run)
            new_run_fl = False
        elif inp == 'q':
            # Go to next run
            ix_run+=1
            new_run_fl = True
        elif inp == 'qq':
            # Quit and save
            next_run_fl = True
            quit_fl = True
        else:
            print('Action not recognized')
            new_run_fl = False

        if ix_run > (len(runs) - 1):
            next_run_fl = True
            pst[row.Index] = True

    if quit_fl:
        break
    
    plt.clf()
        # Ask to re-fit



        # If re-fit, replot


        # Ask if want to apply current re-fit params to open emacs file



        # Ask if want to reset re-fit params (read from init)



        # Ask if want to add multi fit
        # If so, initiate multi fit
        # Display



        # Ask if want to go back (within object)


        # Calculate all sizes


        # Ask which fit to use? (JWST or HST, single or multi, etc)


        # Save

pst.to_pickle(pdata_path + '/galfit/iterfit_status.pkl')
pco.to_pickle(pdata_path + '/galfit/iterfit_comment.pkl')
pcl.to_pickle(pdata_path + '/galfit/iterfit_complab.pkl')
plt.close()

# Object number: 
## 0) sky                    #  object type
## 1) 0.2      1          #  sky background at center of fitting region [ADUs]
## 2) 0.0000      0          #  dsky/dx (sky gradient in x)
## 3) 0.0000      0          #  dsky/dy (sky gradient in y)
## Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 
