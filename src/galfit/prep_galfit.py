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

def galfit(fin):
    cmd = '/home/tonypahl/galfit/galfit {}'.format(fin)
    sp = subprocess.Popen(cmd, shell=True)
    sp.wait()
    return cmd + '\n'

def invert_sqrt_fits(arr_in,arr_out):

    with fits.open(arr_in) as hdu:
        data = hdu[0].data
        data2 = np.sqrt(1 / data)
        hdu[0].data = data2

        hdu.writeto(arr_out,overwrite=True)

# Retrieve CEERS objects to fit
# Use VDW fits as initial guess
#This pickled Pandas data frame is then assigned to pcrs
pcrs = pd.read_pickle('{}/jc_xi_ion_prop.pkl'.format(data_path))

pcrs = pcrs.loc[pcrs.index.str.contains('JADES')]


if len(ix_tofit) > 0:
    pcrs = pcrs.loc[ix_tofit]

for row in pcrs.itertuples():
    for run in run_strs:
        fi = fi_dir[run]
        field = 'GOODS-S'

        # Ignore VDW fit, if no match
        if 'vdw' in run and (pd.isnull(row.ID_3DHST)):
            continue
        
        # Check field membership
        if pd.isnull(field):
            continue
        pxsc = pxsc_dir[fi][field]
        
        # Create directories and link stamps, psf
        id_dir = pdata_path + '/galfit/{}'.format(row.Index)
        if not os.path.isdir(id_dir):
            os.mkdir(id_dir)
        run_dir = id_dir + '/{}'.format(run)
        # Ignore multi runs that haven't been initialized
        if (('multi' in run) and
            (not os.path.isfile(run_dir + '/{}.{}.fit.fits'.format(row.Index, fi)))):
            continue
        if not os.path.isdir(run_dir):
            os.mkdir(run_dir)
        stamp_orig_dir = pdata_path + '/stamps/{}/{}'.format(row.Index, fi)
        stamp_dir = run_dir + '/stamps'
        if not os.path.islink(stamp_dir):
            os.symlink(stamp_orig_dir, stamp_dir, target_is_directory=True)
        #psf_orig_file = pdata_path + '/psf/{}_{:03.0f}mas_psf.fits'.format(fi, pxsc*1000)
        psf_orig_file = pdata_path + '/psf/{}_{}_psf.fits'.format(fi, field)
        psf_file_wpath = run_dir + '/{}.{}.psf.fits'.format(row.Index, run)
        if os.path.islink(psf_file_wpath):
            os.unlink(psf_file_wpath)
        os.symlink(psf_orig_file, psf_file_wpath)

        # Don't initilize multi runs further, as iterfit takes care of that
        if 'multi' in run:
            continue

        # Create galfit config file
        with open('readable_in.galfit','r') as f:
            config_lines = f.readlines()

        # Configure images
        sci_file = './stamps/{}_{}_sci.2as.fits'.format(row.Index, fi)
        out_file = './{}.{}.fit.fits'.format(row.Index, run)
        bpm_file = './stamps/{}_{}_bpm.2as.fits'.format(row.Index, fi)
        psf_file = './{}'.format(os.path.basename(psf_file_wpath))
        sig_file = './stamps/{}_{}_sig.2as.fits'.format(row.Index, fi) #'none'
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
        if not os.path.isfile('{}/{}'.format(run_dir, sci_file)):
            print(row.Index)
            continue

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
        config_lines = [subs.replace('!HOW!',str(how_dir[run])) for subs in config_lines]
        how_arr = how_param_dir[run]
        

        # Sersic profiles
        # Define x, y positions, and other parameters from vdw, if available
        if False: #pd.isnull(row.ID_3DHST):
            wcs = WCS('{}/{}'.format(run_dir, sci_file))
            x_image, y_image = wcs.world_to_pixel(SkyCoord(row.ra_vdw, row.dec_vdw, unit='deg'))
            n = row.n_vdw
            r_e = row.re_vdw / pxsc # r_vdw is in arcsec
            pa = row.pa_vdw
            mag = row.mag_vdw
            q = row.q_vdw
            # If VDW has optimized, use that x, y position
            if 'vdw' not in run:
                vdw_fit = id_dir + '/F160W_vdw/{}.F160W_vdw.fit.fits'.format(row.Index)
                if os.path.isfile(vdw_fit):
                    with fits.open(vdw_fit) as hdu:
                        model_header = hdu[2].header
                        x_image = float(model_header['1_XC'].split(' ')[0].replace('*','').replace('[','').replace(']',''))
                        y_image = float(model_header['1_YC'].split(' ')[0].replace('*','').replace('[','').replace(']',''))
                    if not (fi == 'F160W'):
                        field_F160W = getattr(row, 'field_F160W')
                        xy_scale = pxsc_dir['F160W'][field_F160W] / pxsc_dir['F115W'][field]
                        x_image= x_image * xy_scale
                        y_image = y_image * xy_scale
        elif False: #pd.isnull(getattr(row, '{}_SEx_ID'.format(fi))): # Otherwise, use SEx catalog, if available
            x_image, y_image = (np.array(stamp_shape) - 1) / 2.
            n = 4
            r_e = getattr(row, 're_{}'.format(fi))
            pa = 0
            mag = getattr(row, 'mag_{}'.format(fi))
            q = 1
        else: # Otherwise, use default settings
            x_image, y_image = (np.array(stamp_shape) - 1) / 2.
            n = 4
            r_e = 3.
            pa = 0
            if fi == 'F160':
                mag = 27.
            else:
                mag = 28.
            q = 1
            
                    

        # Constraints (if any)
        const_str = constraints_dir[run]
        with open('{}/galfit.constraints'.format(run_dir), 'w') as f:
            f.write(const_str)

        sersic = sersic_str(1, x_image, y_image, n=n, r_e=r_e, pa=pa, mag=mag, q=q, how=how_arr)
        single_param = copy(config_lines)
        single_param.extend(sersic)
        single_param.extend(end_str())
        with open(run_dir + '/{}.{}.init.galfit'.format(row.Index, run), 'w') as f:
            f.writelines(single_param)

        # If it doesn't exist, create "run" file
        run_file = run_dir + '/{}.{}.run.galfit'.format(row.Index, run)
        if not os.path.isfile(run_file):
            with open(run_file, 'w') as f:
                f.writelines(single_param)
        
        
        
        
        
        
        
