
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

from vars import *

def galfit(fin):
    cmd = '/home/tonypahl/galfit/galfit {}'.format(fin)
    sp = subprocess.Popen(cmd, shell=True)
    sp.wait()
    return cmd + '\n'


def refit_galfit(obj, run):
    init_dir = os.getcwd() + '/'
    
    id_dir = pdata_path + '/galfit/{}'.format(obj)
    run_dir = id_dir + '/{}'.format(run)

    
    os.chdir(run_dir)
    arg = '{}.{}.run.galfit'.format(obj, run)
    res = galfit(arg)
    print(res)

    os.chdir(init_dir)
