# Written by Anthony Pahl, 2/16/23

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

from astropy.wcs.utils import skycoord_to_pixel
import sewpy

from vars import *

def invert_sqrt_fits(arr_in,arr_out):

    with fits.open(arr_in) as hdu:
        data = hdu[0].data
        data2 = np.sqrt(1 / data)
        hdu[0].data = data2

        hdu.writeto(arr_out,overwrite=True)

pcrs = pd.read_pickle('{}/cat/ceers_radec.pkl'.format(res_path))

# Read in all science images
sci_img_dir = copy.deepcopy(sci_dir)
sci_header_dir = copy.deepcopy(sci_dir)
sci_wcs_dir = copy.deepcopy(sci_dir)
for fi in sci_dir:
    for field in sci_dir[fi]:
        with fits.open(sci_dir[fi][field]) as hdu:
            img = hdu[0].data
            wcs = WCS(hdu[0].header)
            header = hdu[0].header
        sci_img_dir[fi][field] = img
        sci_header_dir[fi][field] = header
        sci_wcs_dir[fi][field] = wcs

# Read in all weight images
wht_img_dir = copy.deepcopy(wht_dir)
wht_header_dir = copy.deepcopy(wht_dir)
wht_wcs_dir = copy.deepcopy(wht_dir)
for fi in wht_dir:
    for field in wht_dir[fi]:
        with fits.open(wht_dir[fi][field]) as hdu:
            img = hdu[0].data
            wcs = WCS(hdu[0].header)
            header = hdu[0].header
        wht_img_dir[fi][field] = img
        wht_header_dir[fi][field] = header
        wht_wcs_dir[fi][field] = wcs

# Read in all segmentation images
seg_img_dir = copy.deepcopy(seg_dir)
seg_header_dir = copy.deepcopy(seg_dir)
seg_wcs_dir = copy.deepcopy(seg_dir)
for field in seg_dir:
    with fits.open(seg_dir[field]) as hdu:
        img = hdu[0].data
        wcs = WCS(hdu[0].header)
        header = hdu[0].header
    seg_img_dir[field] = img
    seg_header_dir[field] = header
    seg_wcs_dir[field] = wcs

# Determine image membership
for row in pcrs.itertuples():
    coord = SkyCoord(row.ra, row.dec, unit=u.deg)
    for fi in filters:
        field_mem = np.empty_like(field_ord, dtype=bool)
        for i, field in enumerate(field_ord):
            if field not in sci_dir[fi]:
                #pcrs.loc[row.Index, 'field_{}'.format(fi)] = np.nan
                field_mem[i] = False
            else:
                if not coord.contained_by(sci_wcs_dir[fi][field]):
                    field_mem[i] = False
                else:
                    x,y = skycoord_to_pixel(coord, sci_wcs_dir[fi][field])
                    if sci_img_dir[fi][field][int(y), int(x)] == 0:
                        field_mem[i] = False
                    else:
                        field_mem[i] = True
        field_oks = [i for (i, v) in zip(field_ord, field_mem) if v]
        if len(field_oks) == 0:
            pcrs.loc[row.Index, 'field_{}'.format(fi)] = np.nan
        else:
            pcrs.loc[row.Index, 'field_{}'.format(fi)] = field_oks[0]
pcrs.loc[:, ['field_F115W', 'field_F160W']].to_pickle('{}/stamps/field_mem.pkl'.format(pdata_path))

# Cutout stamps
for row in pcrs.itertuples():
    coord = SkyCoord(row.ra, row.dec, unit=u.deg)
    workdir = '{}/stamps/{}'.format(pdata_path, row.Index)
    if not os.path.isdir(workdir):
        os.mkdir(workdir)

    for fi in filters:
        field = pcrs.loc[row.Index, 'field_{}'.format(fi)]
        if pd.isnull(field):
            continue
        
        workdir = '{}/stamps/{}/{}'.format(pdata_path, row.Index, fi)
        if not os.path.isdir(workdir):
            os.mkdir(workdir)
        
        # Define size of all stamps
        st_size = 5.0 * u.arcsec #arcsec
        px_size = st_size  / pxsc_dir[fi][field]
        
        # Science
        sci_stamp_file = '{}/{}_{}_sci.5as.fits'.format(workdir, row.Index, fi)
        cutout_sci = Cutout2D(
            sci_img_dir[fi][field],
            coord,
            st_size,
            wcs=sci_wcs_dir[fi][field]
        )
        hdu_sci = fits.PrimaryHDU()
        hdu_sci.data = cutout_sci.data
        hdu_sci.header = sci_header_dir[fi][field]
        hdu_sci.header.update(cutout_sci.wcs.to_header())
        if 'EXPTIME' in hdu_sci.header:
            hdu_sci.header.update(EXPTIME=1)
        hdu_sci.writeto(sci_stamp_file, overwrite=True)

        # Weight
        wht_stamp_file = '{}/{}_{}_wht.5as.fits'.format(workdir, row.Index, fi)
        cutout_wht = Cutout2D(
            wht_img_dir[fi][field],
            cutout_sci.position_original,
            cutout_sci.shape,    
            wcs=wht_wcs_dir[fi][field]
        )
        hdu_wht = fits.PrimaryHDU()
        hdu_wht.data = cutout_wht.data
        hdu_wht.header = wht_header_dir[fi][field]
        hdu_wht.header.update(cutout_wht.wcs.to_header())
        hdu_wht.writeto(wht_stamp_file, overwrite=True)

        # Sigma
        sig_stamp_file = '{}/{}_{}_sig.5as.fits'.format(workdir, row.Index, fi)
        invert_sqrt_fits(wht_stamp_file, sig_stamp_file)

        # Segmentation map
        # Check if object is imaged in F115W, if not, use custom generated segmentation map
        if not (fi == 'F115W') and pd.isnull(getattr(row, 'field_F115W')):
            segfield = 'AEGIS'
        else:
            segfield = field            
        seg_stamp_file = '{}/{}_{}_seg.5as.fits'.format(workdir, row.Index, fi)
        # Create segmap by checking where each pixel lies in seg image (works for different
        # pixel scales)
        x2D, y2D = np.meshgrid(np.arange(cutout_sci.shape[0]), np.arange(cutout_sci.shape[1]))
        xys = np.column_stack((y2D.ravel(), x2D.ravel()))
        coords = cutout_sci.wcs.pixel_to_world(xys[:,0], xys[:,1])
        ixs = seg_wcs_dir[segfield].world_to_array_index(coords)
        seg_data = seg_img_dir[segfield][ixs[0], ixs[1]].reshape(cutout_sci.data.shape).T
        hdu_seg = fits.PrimaryHDU()
        hdu_seg.data = seg_data
        hdu_seg.header = seg_header_dir[segfield]
        hdu_seg.header.update(cutout_sci.wcs.to_header())
        hdu_seg.writeto(seg_stamp_file, overwrite=True)

        # Bad pixel map
        bpm_stamp_file = '{}/{}_{}_bpm.5as.fits'.format(workdir, row.Index, fi)
        # Create bpm by replacing central value by 0
        cen_val = seg_data[int(seg_data.shape[0]/2), int(seg_data.shape[0]/2)]
        bpm_data = copy.copy(seg_data)
        bpm_data[bpm_data == cen_val] = 0
        #bpm_data[bpm_data > 0] = 1
        hdu_bpm = fits.PrimaryHDU()
        hdu_bpm.data = bpm_data
        hdu_bpm.header = seg_header_dir[segfield]
        hdu_bpm.header.update(cutout_sci.wcs.to_header())
        hdu_bpm.writeto(bpm_stamp_file, overwrite=True)
            

