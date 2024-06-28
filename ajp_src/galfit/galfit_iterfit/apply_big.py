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

def apply_big(obj, fi, row, run):
    
    id_dir = pdata_path + '/galfit/{}'.format(obj)
    run_dir = id_dir + '/{}'.format(run)

    # Open currently active runfile
    run_file = run_dir + '/{}.{}.run.galfit'.format(row.Index, run)

    with open(run_file) as f:
        config_lines = f.readlines()

    sci_file = './stamps/{}_{}_sci.fits'.format(row.Index, fi)
    bpm_file = './stamps/{}_{}_bpm.fits'.format(row.Index, fi)
    sig_file = './stamps/{}_{}_sig.fits'.format(row.Index, fi) #'none'

    nsci_file = './stamps/{}_{}_sci.5as.fits'.format(row.Index, fi)
    nbpm_file = './stamps/{}_{}_bpm.5as.fits'.format(row.Index, fi)
    nsig_file = './stamps/{}_{}_sig.5as.fits'.format(row.Index, fi) #'none'
    
    
    config_lines = [subs.replace(sci_file,nsci_file) for subs in config_lines]
    config_lines = [subs.replace(bpm_file,nbpm_file) for subs in config_lines]
    config_lines = [subs.replace(sig_file,nsig_file) for subs in config_lines]

    # Define central 200px
    with fits.open('{}/{}'.format(run_dir, nsci_file)) as hdu:
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

    xylim_line_ix = [i for i,x in enumerate(config_lines) if x[:2]=='H)'][0]
    nxylim_line = 'H)\t{}\t{}\t{}\t{}\t# Image region to fit (xmin xmax ymin ymax)\n'.format(
        llims[0], ulims[0], llims[1], ulims[1])
    config_lines[xylim_line_ix] = nxylim_line

    # Fix positions
    xy_line_ix = [i for i,x in enumerate(config_lines) if x[:3]==' 1)']
    if len(xy_line_ix) == 1:
        for ix in xy_line_ix:
            slines = config_lines[ix].split()
            x1 = slines[1]
            x2 = slines[2]
            slines[1] = str(float(x1) * 2)
            slines[2] = str(float(x2) * 2)
            nslines = ' '.join(slines)
            config_lines[ix] = ' ' + nslines + '\n'
        
        
    

    with open(run_file, 'w') as f:
        f.writelines(config_lines)
