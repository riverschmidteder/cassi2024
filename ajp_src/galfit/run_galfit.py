
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

def galfit(fin):
    cmd = '/home/tonypahl/galfit/galfit {}'.format(fin)
    sp = subprocess.Popen(cmd, shell=True)
    sp.wait()
    return cmd + '\n'

# Retrieve CEERS objects to fit
pcrs = pd.read_pickle(res_path + '/cat/ceers_vdw.pkl')
pphot = pd.read_pickle(res_path + '/cat/ceers_sex_phot.pkl')
pcrs = pphot.merge(pcrs.iloc[:,4:], on='ID', how='left')

if len(ix_tofit) > 0:
    pcrs = pcrs.loc[ix_tofit]

init_dir = os.getcwd() + '/'

for row in pcrs.itertuples():
    for run in run_strs:
        id_dir = pdata_path + '/galfit/{}'.format(row.Index)
        run_dir = id_dir + '/{}'.format(run)

        # Check field membership
        fi = fi_dir[run]
        field = getattr(row, 'field_{}'.format(fi))
        if pd.isnull(field):
            continue
        if not os.path.isdir(run_dir):
            continue

        # Run
        os.chdir(run_dir)
        arg = '{}.{}.run.galfit'.format(row.Index, run)
        res = galfit(arg)
        print(res)


os.chdir(init_dir)

        
        
        
        
        
        
