import os, subprocess, sys, shutil
import glob
import pandas as pd
from pandas import Series, DataFrame
from astropy.table import Table
import numpy as np
from copy import copy
import multiprocessing as mp
import subprocess

from astropy.io import fits,ascii
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.nddata import Cutout2D
import astropy.units as u

from vars import *

def sersic_str(num, x, y, mag=20, r_e=5, n=4, q=0.6, pa=60, how=[1, 1, 1, 1, 1, 1, 1]):
    hx, hy, hmag, hr_e, hn, hq, hpa = how
    out_strs = [
        '# Object number: {}\n'.format(int(num)),
        ' 0) sersic                 #  object type\n',
        ' 1) {}  {}  {} {}  #  position x, y\n'.format(x,y,hx,hy),
        ' 3) {}     {}          #  Integrated magnitude\t\n'.format(mag,hmag),
        ' 4) {}      {}          #  R_e (half-light radius)   [pix]\n'.format(r_e,hr_e),
        ' 5) {}      {}          #  Sersic index n (de Vaucouleurs n=4) \n'.format(n,hn),
        ' 6) 0.0000      0          #     ----- \n',
        ' 7) 0.0000      0          #     ----- \n',
        ' 8) 0.0000      0          #     ----- \n',
        ' 9) {}      {}          #  axis ratio (b/a)  \n'.format(q,hq),
        '10) {}    {}          #  position angle (PA) [deg: Up=0, Left=90]\n'.format(pa,hpa),
        " Z) 0                      #  output option (0 = resid., 1 = Don't subtract) \n",
        '\n',
        ]
    return out_strs

def sky_str(num):
    out_strs = [
        '# Object number: {}\n'.format(num),
        ' 0) sky                    #  object type\n',
        ' 1) 1.3920      1          #  sky background at center of fitting region [ADUs]\n',
        ' 2) 0.0000      0          #  dsky/dx (sky gradient in x)\n',
        ' 3) 0.0000      0          #  dsky/dy (sky gradient in y)\n',
        " Z) 0                      #  output option (0 = resid., 1 = Don't subtract) \n",
        '\n',
        ]
    return out_strs

def end_str():
    return [ '================================================================================\n',
 '\n']


def init_mult_fit(obj, fi, row):
    run = '{}_multi'.format(fi)
    field = getattr(row, 'field_{}'.format(fi))
    pxsc = pxsc_dir[fi][field]
    
    # Create directories and link stamps, psf
    id_dir = pdata_path + '/galfit/{}'.format(obj)
    run_dir = id_dir + '/{}'.format(run)
    if not os.path.isdir(run_dir):
        os.mkdir(run_dir)
    stamp_orig_dir = pdata_path + '/stamps/{}/{}'.format(row.Index, fi)
    stamp_dir = run_dir + '/stamps'
    if not os.path.islink(stamp_dir):
        os.symlink(stamp_orig_dir, stamp_dir, target_is_directory=True)
    psf_orig_file = pdata_path + '/psf/{}_{}_psf.fits'.format(fi, field)
    psf_file_wpath = run_dir + '/{}.{}.psf.fits'.format(row.Index, run)
    if not os.path.islink(psf_file_wpath):
        os.symlink(psf_orig_file, psf_file_wpath)
    
        
    # Create galfit config file
    with open('readable_in.galfit','r') as f:
        config_lines = f.readlines()

    # Configure images
    sci_file = './stamps/{}_{}_sci.fits'.format(row.Index, fi)
    out_file = './{}.{}.fit.fits'.format(row.Index, run)
    bpm_file = './stamps/{}_{}_bpm.fits'.format(row.Index, fi)
    psf_file = './{}'.format(os.path.basename(psf_file_wpath))
    sig_file = './stamps/{}_{}_sig.fits'.format(row.Index, fi) #'none'
    config_lines = [subs.replace('!SCIENCE!',sci_file) for subs in config_lines]
    config_lines = [subs.replace('!OUTPUT!',out_file) for subs in config_lines]
    config_lines = [subs.replace('!PSF!',psf_file) for subs in config_lines]
    config_lines = [subs.replace('!BPM!',bpm_file) for subs in config_lines]
    config_lines = [subs.replace('!SIGMA!',sig_file) for subs in config_lines]

    # Zerpoint
    mag_zpt = zeropoint_dir[fi][field]
    config_lines = [subs.replace('!ZEROPT!',str(mag_zpt)) for subs in config_lines] #temporary

    # Plate scale
    config_lines = [subs.replace('!PS!',str(pxsc)) for subs in config_lines] #temporary


    # Define central 200px
    with fits.open('{}/{}'.format(run_dir, sci_file)) as hdu:
        stamp_shape = hdu[0].data.shape
    llims = np.empty(2,dtype='object')
    ulims = np.empty(2,dtype='object')
    for i in range(2):
        imin = stamp_shape[i]/2. - 100
        imax = stamp_shape[i]/2. + 100

        if imin < 1:
            llims[i] = '1'
            ulims[i] = str(stamp_shape[i]-2)
        else:
            llims[i] = str(int(np.ceil(imin)))
            ulims[i] = str(int(np.floor(imax)))
    config_lines = [subs.replace('!XMIN!',llims[0]) for subs in config_lines]
    config_lines = [subs.replace('!XMAX!',ulims[0]) for subs in config_lines]
    config_lines = [subs.replace('!YMIN!',llims[1]) for subs in config_lines]
    config_lines = [subs.replace('!YMAX!',ulims[1]) for subs in config_lines]

    # Choose to either optimize or output model
    config_lines = [subs.replace('!HOW!',str(0)) for subs in config_lines]
    #how_arr = how_param_dir[run]


    # Sersic profiles
    # Define x, y positions, and other parameters from best-fit single
    fit_file_single = '{}/galfit/{}/{}/{}.{}.fit.fits'.format(
        pdata_path,
        row.Index,
        run.replace('multi','single'),
        row.Index,
        run.replace('multi','single'))
    with fits.open(fit_file_single) as hdu:
        model_header = hdu[2].header
    x_image = float(model_header['1_XC'].split(' ')[0].replace('*','').replace('[','').replace(']',''))
    y_image = float(model_header['1_YC'].split(' ')[0].replace('*','').replace('[','').replace(']',''))
    mag = float(model_header['1_MAG'].split(' ')[0].replace('*','').replace('[','').replace(']',''))
    r_e = float(model_header['1_RE'].split(' ')[0].replace('*','').replace('[','').replace(']',''))
    n = float(model_header['1_N'].split(' ')[0].replace('*','').replace('[','').replace(']',''))
    pa = float(model_header['1_PA'].split(' ')[0].replace('*','').replace('[','').replace(']',''))
    q = float(model_header['1_AR'].split(' ')[0].replace('*','').replace('[','').replace(']',''))

    # Constraints (if any)
    const_str = constraints_dir[run]
    with open('{}/galfit.constraints'.format(run_dir), 'w') as f:
        f.write(const_str)

    sersic = sersic_str(1, x_image, y_image, n=n, r_e=r_e, pa=pa, mag=mag, q=q)
    sersic2 = sersic_str(2, x_image, y_image, n=n, r_e=r_e, pa=pa, mag=mag, q=q)
    single_param = copy(config_lines)
    single_param.extend(sersic)
    single_param.extend(sersic2)
    single_param.extend(end_str())
    
    with open(run_dir + '/{}.{}.init.galfit'.format(row.Index, run), 'w') as f:
        f.writelines(single_param)

    # If it doesn't exist, create "run" file
    run_file = run_dir + '/{}.{}.run.galfit'.format(row.Index, run)
    if not os.path.isfile(run_file):
        with open(run_file, 'w') as f:
            f.writelines(single_param)
    else:
        run = np.nan

    return run

def remove_mult_fit(obj, fi, row):
    run = '{}_multi'.format(fi)
    field = getattr(row, 'field_{}'.format(fi))
    pxsc = pxsc_dir[fi][field]
    
    # Create directories and link stamps, psf
    id_dir = pdata_path + '/galfit/{}'.format(obj)
    run_dir = id_dir + '/{}'.format(run)
    if not os.path.isdir(run_dir):
        return np.nan
    else:
        shutil.rmtree(run_dir)
        return run
        
        
