# Written by Anthony Pahl, 2/16/23

import os
import subprocess
import sys
import shutil
import glob
import pandas as pd
from pandas import Series, DataFrame
from astropy.table import Table
import numpy as np
import copy

from astropy.io import fits, ascii
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

from vars import *

pcrs = pd.read_pickle('{}/cat/ceers_radec.pkl'.format(res_path))
pfield_mem = pd.read_pickle('{}/stamps/field_mem.pkl'.format(pdata_path))
pcrs = pcrs.merge(pfield_mem, left_index=True, right_index=True)

out_file = '{}/stamps/stamps.pdf'.format(pdata_path)

with PdfPages(out_file) as pdf:
    for row in pcrs.itertuples():
        fis = copy.copy(filters)
        for i, fi in enumerate(filters):
            if pd.isnull(pcrs.loc[row.Index, 'field_{}'.format(fi)]):
                del fis[i]

        num_x = 4
        num_y = len(fis)
        fig = plt.figure(figsize=(4 * num_x + 0.5, 4 * num_y + 1))
        gs = gridspec.GridSpec(num_y + 1, num_x + 1,
                               height_ratios=[1] + [4] * num_y,
                               width_ratios=[0.5] + [4] * num_x)

        # Info
        ax_info = fig.add_subplot(gs[0, 1:])
        ax_info.axis('off')
        ax_info.text(
            0.0, 0.9,
            '{}\nz={:.4f}'.format(row.Index, row.z),
            ha='left', va='top',
            fontsize=20,
            transform=ax_info.transAxes
        )

        for i, fi in enumerate(fis):
            row_pos = 1 + i

            # Science
            col_pos = 1
            ax = fig.add_subplot(gs[row_pos, col_pos])
            sci_stamp_file = '{}/stamps/{}/{}/{}_{}_sci.fits'.format(
                pdata_path,
                row.Index,
                fi,
                row.Index, fi)
            with fits.open(sci_stamp_file) as hdu:
                sci_data = hdu[0].data
            ax.imshow(sci_data, cmap='gist_yarg', origin='lower')
            ax.text(
                -0.1, 0.5, fi,
                size=18,
                ha='center', va='bottom',
                rotation=90., rotation_mode='anchor',
                transform=ax.transAxes)
            if row_pos == 1:
                ax.text(
                    0.5, 1.1, 'Science',
                    size=18,
                    ha='center', va='bottom',
                    transform=ax.transAxes)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

            # Weight
            col_pos = 4
            ax = fig.add_subplot(gs[row_pos, col_pos])
            wht_stamp_file = '{}/stamps/{}/{}/{}_{}_wht.fits'.format(
                pdata_path,
                row.Index,
                fi,
                row.Index, fi)
            with fits.open(wht_stamp_file) as hdu:
                wht_data = hdu[0].data
            ax.imshow(wht_data, cmap='gist_yarg', origin='lower')
            if row_pos == 1:
                ax.text(
                    0.5, 1.1, 'Weight',
                    size=18,
                    ha='center', va='bottom',
                    transform=ax.transAxes)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

            # Segmentation map
            col_pos = 2
            ax = fig.add_subplot(gs[row_pos, col_pos])
            seg_stamp_file = '{}/stamps/{}/{}/{}_{}_seg.fits'.format(
                pdata_path,
                row.Index,
                fi,
                row.Index, fi)
            with fits.open(seg_stamp_file) as hdu:
                seg_data = hdu[0].data
            seg_values = set(seg_data.flatten())
            seg_newvals = np.arange(len(seg_values))
            for i, val in enumerate(seg_values):
                seg_data[seg_data == val] = seg_newvals[i]
            imseg = ax.imshow(
                seg_data,
                cmap='gnuplot',
                origin='lower',
                interpolation='none')
            clim_seg = imseg.properties()['clim']
            # plot central dot
            circle1 = plt.Circle((0.5, 0.5),
                                 0.005, color='w',
                                 transform=ax.transAxes)
            circle2 = plt.Circle((0.5, 0.5),
                                 0.01, color='b',
                                 transform=ax.transAxes)
            ax.add_patch(circle2)
            ax.add_patch(circle1)
            if row_pos == 1:
                ax.text(
                    0.5, 1.1, 'Seg',
                    size=18,
                    ha='center', va='bottom',
                    transform=ax.transAxes)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

            # Bad pixel map
            col_pos = 3
            ax = fig.add_subplot(gs[row_pos, col_pos])
            bpm_stamp_file = '{}/stamps/{}/{}/{}_{}_bpm.fits'.format(
                pdata_path,
                row.Index,
                fi,
                row.Index, fi)
            with fits.open(bpm_stamp_file) as hdu:
                bpm_data = hdu[0].data
            for i, val in enumerate(seg_values):
                bpm_data[bpm_data == val] = seg_newvals[i]
            ax.imshow(bpm_data, cmap=imseg.get_cmap(), clim=clim_seg,
                      origin='lower', interpolation='none')
            if row_pos == 1:
                ax.text(
                    0.5, 1.1, 'Bpm',
                    size=18,
                    ha='center', va='bottom',
                    transform=ax.transAxes)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

        gs.tight_layout(fig, h_pad=0.1, w_pad=0.1)
        pdf.savefig()
        plt.close()
