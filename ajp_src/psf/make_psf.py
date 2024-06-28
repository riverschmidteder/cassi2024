import os, subprocess, sys
import glob
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from astropy.visualization import simple_norm

from astropy.table import Table

from astropy.io import fits,ascii
from stwcs import updatewcs
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

from astropy.nddata import NDData
from photutils.psf import extract_stars
from photutils import EPSFBuilder
from photutils import centroids

from astropy.modeling import models, fitting

from vars import *


def write_region(df,idxs,out, color, radius):
    with open(out,'w') as f:
        f.write('''# Region file format: DS9 version 4.1
        global color={} dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
        fk5
        '''.format(color))
        for i,r in df.iterrows():
            #if r.CLASS_STAR > thresh:
            f.write('circle({},{},{}) # width=3 text={{{}}} \n'.format(r.RA,r.DEC,radius,idxs[i]))


def plot_stars(stars,idx,outfile=None,show=True):
    nrows = 3
    ncols = 5
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 7),
                           squeeze=True)
    ax = ax.ravel()
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    for i in range(stars.n_stars):
        norm = simple_norm(stars[i], 'log', percent=99.)
        ax[i].imshow(stars[i], norm=norm, origin='lower', cmap='viridis')

        
        ax[i].text(0.05, 0.95, idx[i] , transform=ax[i].transAxes,
                verticalalignment='top', bbox=props)

        #fit 2D gaussian for FWHM
        data = stars[i].data
        p_init = models.Gaussian2D(np.max(data),data.shape[1],data.shape[0],( 4 / 2.355), (4 / 2.355))
        fit_p = fitting.LevMarLSQFitter()

        yi, xi = np.indices(data.shape)
        g = fit_p(p_init, xi, yi, data)
        #print(p_init)
        model_data = g(xi, yi)
        fwhm = (g.x_fwhm + g.y_fwhm) / 2.

        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        #ax[i].text(0.05, 0.85,  'FWHM: {:.2f}'.format(fwhm), transform=ax[i].transAxes,
        #            verticalalignment='top', bbox=props)

    fig.tight_layout()

    if outfile:
        fig.savefig(outfile)

    if show:
        plt.show()

    plt.close(fig)

def plot_psf(epsf,nfilt,outfile=None,show=True):

    norm = simple_norm(epsf.data, 'log', percent=99.)
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111)
    cax = ax.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
    fig.colorbar(cax,ax=ax)

    #fit 2D Gaussian
    p_init = models.Gaussian2D(np.max(epsf.data),epsf.data.shape[1],epsf.data.shape[0],( 6 / 2.355), (6 / 2.355))
    fit_p = fitting.LevMarLSQFitter()

    yi, xi = np.indices(epsf.data.shape)
    g = fit_p(p_init, xi, yi, epsf.data)
    #print(p_init)
    model_data = g(xi, yi)
    fwhm = (g.x_fwhm + g.y_fwhm) / 2.
    
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    #ax.text(0.05, 0.95,  'FWHM: {:.2f}'.format(fwhm), transform=ax.transAxes,
    #            verticalalignment='top', bbox=props)

    ax.set_title(nfilt)
    
    if outfile:
        fig.savefig(outfile)

    if show:
        plt.show()

    plt.close(fig)

for field in fields:

    for fi in filters:

        if fi == 'F115W' and field == 'AEGIS':
            continue

        starlist = starlist_dir[fi][field]

        outstr = '{}/{}_{}'.format(pdata_path, fi, field)

        image = sci_dir[fi][field]
        
        hdu = fits.open(image)
        im = hdu[0].data
        wcs = WCS(hdu[0].header)

        nddata = NDData(data=im)

        pstar = ascii.read(starlist,names=['x','y']).to_pandas()

        out_ix = []
        for i,r in pstar.iterrows():
            if np.isnan(im[int(r.y),int(r.x)]):
                out_ix.append(i)
        pstar = pstar.reindex(index = pstar.index.drop(out_ix))
        atstar = Table.from_pandas(pstar)
                
            
        

        stars = extract_stars(nddata,atstar,69)
        pstarobj = pd.DataFrame(data={'star_obj': np.array(stars.all_stars), 'flux': stars.flux}, dtype='object')
        ##remove galaxies / bad stars
        #if dirx == 'Q0933':
        #    if fi is not 'F606W':
        #        keep_idx = tokeep_code[dirx][fi]
        #        pstarobj = pstarobj.loc[keep_idx]
        #    else:
        #        drop_idx = todrop_code[dirx][fi]
        #        pstarobj = pstarobj.reindex(index = pstarobj.index.drop(drop_idx))
        #else:
        #drop_idx = todrop_code[dirx][fi]
        #pstarobj = pstarobj.reindex(index = pstarobj.index.drop(drop_idx))
                
        
        #choose 15 brightest stars
        pstarobj = pstarobj.sort_values(by='flux',ascending=False)
        bright_idx = pstarobj.iloc[:15].index.values
        #if fi == 'F115W':
        #    raise Exception()


        stars = extract_stars(nddata,atstar[bright_idx],69)



        plot_stars(stars,bright_idx,outfile=outstr+'_stars.png',show=True)

        epsf_builder = EPSFBuilder(oversampling=1,maxiters = 10, recentering_func=centroids.centroid_com)
        epsf,fitted_stars = epsf_builder(stars)

        
        #write stars to region file
        outn_reg = outstr + '_psf_stars.reg'
        preg = pd.DataFrame(fitted_stars.center_flat, columns=['X_IMAGE','Y_IMAGE'])
        preg_radec = pd.DataFrame().reindex_like(preg)
        preg_radec.columns = ['RA', 'DEC']
        for i,r in preg.iterrows():
            coord = wcs.pixel_to_world(r.X_IMAGE, r.Y_IMAGE)
            preg_radec.loc[i, 'RA'] = coord.ra.value
            preg_radec.loc[i, 'DEC'] = coord.dec.value

       ## Write out starlist for F160W
       #if fi == 'F115W':
       #    starlist_F160W = starlist_dir['F160W'][field]
       #
       #    image = sci_dir['F160W'][field]
       #
       #    hdu = fits.open(image)
       #    wcs = WCS(hdu[0].header)
       #    preg_f160w = pd.DataFrame().reindex_like(preg)
       #    for i,r in preg_radec.iterrows():
       #        coord = SkyCoord(r.RA, r.DEC, unit=u.deg)
       #        x, y = wcs.world_to_pixel(coord)
       #        preg_f160w.loc[i, 'X_IMAGE'] = x
       #        preg_f160w.loc[i, 'Y_IMAGE'] = y
       #    preg_f160w.to_csv(starlist_F160W, sep=' ', columns=['X_IMAGE','Y_IMAGE'], header=False, index=False)


            ## Do AEGIS
        #if fi == 'F160W':
        #    starlist_AEGIS = starlist_dir['F160W']['AEGIS']
        #    image = sci_dir['F160W']['AEGIS']
        #    hdu = fits.open(image)
        #    wcs = WCS(hdu[0].header)
        #    preg_f160w = pd.DataFrame().reindex_like(preg)
        #    for i,r in preg_radec.iterrows():
        #        coord = SkyCoord(r.RA, r.DEC, unit=u.deg)
        #        x, y = wcs.world_to_pixel(coord)
        #        preg_f160w.loc[i, 'X_IMAGE'] = x
        #        preg_f160w.loc[i, 'Y_IMAGE'] = y
        #    preg_f160w.to_csv(starlist_AEGIS, mode='a', sep=' ', columns=['X_IMAGE','Y_IMAGE'], header=False, index=False)
            

            
            
        
        write_region(preg_radec,bright_idx,outn_reg,'green','1.0"')

        #write stars to xy
        preg['ix'] = bright_idx
        preg.to_csv(outstr + '_psf_stars.xy',sep = ' ',columns=['X_IMAGE','Y_IMAGE','ix'],header=False,index=False)

        plot_psf(epsf,fi,outfile=outstr+'_psf.png',show=True)


        hdu1 = fits.PrimaryHDU(epsf.data)
        hdunew = fits.HDUList([hdu1])
        hdunew.writeto(outstr+'_psf.fits',overwrite=True)
