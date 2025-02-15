
import os, subprocess, sys, shutil
import glob
import pandas as pd
from pandas import Series, DataFrame
from astropy.table import Table
import numpy as np
import copy

from astropy.io import fits,ascii
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.nddata import Cutout2D
import astropy.units as u
import regions

import numpy.ma as ma

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.backends.backend_pdf import PdfPages
from astropy.visualization import SqrtStretch
from astropy.visualization import ZScaleInterval
from astropy.visualization import make_lupton_rgb

from matplotlib.colors import LogNorm
import matplotlib as mpl

def plot_obj(obj, runs, rinfo, fig):

    
    sys.path.append("..")
    from vars import fi_dir, pxsc_dir, pdata_path, run_strs
    
    num_x = 5
    num_y = len(run_strs)
    #fig = plt.figure(figsize = (4*num_x+0.5, 4*num_y+1))
    gs = gridspec.GridSpec(num_y+1, num_x+1,
                           height_ratios=[1]+[4]*num_y,
                           width_ratios=[0.5]+[4]*num_x,
                           )

    # Info
    ax_info = fig.add_subplot(gs[0,1:])
    ax_info.axis('off')
    ax_info.text(
        0.0, 0.9,
        '{}\nz={:.4f}'.format(rinfo.Index, rinfo.z),
        ha='left', va='top',
        fontsize=20,
        transform=ax_info.transAxes
        )

    for i,run in enumerate(runs):
        row_pos = 1 + i
        fi = fi_dir[run]
        field = getattr(rinfo, 'field_{}'.format(fi))
        pxsc = pxsc_dir[fi][field]


        # Read in Model and Residual information
        fit_file = '{}/galfit/{}/{}/{}.{}.fit.fits'.format(
            pdata_path,
            obj,
            run,
            obj,
            run)
        with fits.open(fit_file) as hdu:
            model_data = hdu[2].data
            model_header = hdu[2].header
            resid_data = hdu[3].data
            model_lower = resid_data.min()
            #model_lower = -1 * (model_data.max() - model_data.min()) * 0.1


        # Science
        col_pos = 2
        ax = fig.add_subplot(gs[row_pos, col_pos])
        sci_stamp_file = '{}/galfit/{}/{}/{}'.format(
            pdata_path,
            obj,
            run,
            model_header['DATAIN'])
        with fits.open(sci_stamp_file) as hdu:
            sci_data = hdu[0].data
        imsci = ax.imshow(
            sci_data,
            cmap='gist_yarg',
            #cmap=imsci.get_cmap(),
            #clim=clim_imsci,
            origin='lower')
        clim_imsci = imsci.properties()['clim']
        clim_imsci = (model_lower, clim_imsci[1])
        imsci = ax.imshow(
            sci_data,
            cmap=imsci.get_cmap(),
            clim=clim_imsci,
            origin='lower')
        if row_pos == 1:
            ax.text(
                0.5, 1.1, 'Science',
                size=18,
                ha='center',va='bottom',
                transform=ax.transAxes)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)


        # Bad pixel map
        col_pos = 1
        ax = fig.add_subplot(gs[row_pos, col_pos])
        bpm_stamp_file = '{}/galfit/{}/{}/{}'.format(
            pdata_path,
            obj,
            run,
            model_header['MASK'])
        if not os.path.isfile(bpm_stamp_file):
            bpm_stamp_file = '{}/galfit/{}/{}/stamps/{}_{}_bpm.fits'.format(
                pdata_path,
                obj,
                run,
                obj, fi)
        with fits.open(bpm_stamp_file) as hdu:
            bpm_data = hdu[0].data
        bpm_values = set(bpm_data.flatten())
        bpm_newvals = np.arange(len(bpm_values))
        for i, val in enumerate(bpm_values):
            bpm_data[bpm_data == val] = bpm_newvals[i]
        ax.imshow(bpm_data, cmap='gnuplot',
                  origin='lower', interpolation='none')
        if row_pos == 1:
            ax.text(
                0.5, 1.1, 'Bpm',
                size=18,
                ha='center',va='bottom',
                transform=ax.transAxes)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        # Run label
        ax.text(
            -0.1, 0.5, run,
            size=18,
            ha='center',va='bottom',
            rotation=90., rotation_mode='anchor',
            transform=ax.transAxes)

        # Model
        col_pos = 3
        ax = fig.add_subplot(gs[row_pos, col_pos])
        ax.imshow(
            model_data,
            cmap=imsci.get_cmap(),
            clim=clim_imsci,
            origin='lower')
        #clim_imsci = imsci.properties()['clim']
        #clim_imsci = (clim_imsci[0] - (clim_imsci[1] - clim_imsci[0])*0.1,
        #              clim_imsci[1])
        if row_pos == 1:
            ax.text(
                0.5, 1.1, 'Model',
                size=18,
                ha='center',va='bottom',
                transform=ax.transAxes)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        # Residual
        col_pos = 4
        ax = fig.add_subplot(gs[row_pos, col_pos])
        ax.imshow(
            resid_data,
            cmap=imsci.get_cmap(),
            clim=clim_imsci,
            origin='lower')
        if row_pos == 1:
            ax.text(
                0.5, 1.1, 'Residual',
                size=18,
                ha='center',va='bottom',
                transform=ax.transAxes)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        # Fit info
        col_pos = 5
        ax = fig.add_subplot(gs[row_pos, col_pos])
        mag = float(model_header['1_MAG'].split(' ')[0].replace('*','').replace('[','').replace(']',''))
        if '1_RE' in model_header:
            r_e = float(model_header['1_RE'].split(' ')[0].replace('*','').replace('[','').replace(']','')) * pxsc * 1000.
        else:
            r_e = np.nan
        if '1_N' in model_header:
            n = float(model_header['1_N'].split(' ')[0].replace('*','').replace('[','').replace(']',''))
        else:
            n = np.nan
        if '1_PA' in model_header:
            pa = float(model_header['1_PA'].split(' ')[0].replace('*','').replace('[','').replace(']',''))
        else:
            pa = np.nan
        if '1_AR' in model_header:
            q = float(model_header['1_AR'].split(' ')[0].replace('*','').replace('[','').replace(']',''))
        else:
            q = np.nan

        if 'label_{}'.format(run) in rinfo._fields:
            label = getattr(rinfo, 'label_{}'.format(run))
        else:
            label = np.nan
        ax.text(0.1, 0.9,
                'mag={:.1f}\nr_e={:.2f}mas\nn={:.1f}\npa={:.1f}\nq={:.2f}\n{}'.format(mag, r_e, n, pa, q,label),
                ha='left', va='top',
                transform=ax.transAxes,
                size=18,
                )
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

    fig.tight_layout()
    #return fig
