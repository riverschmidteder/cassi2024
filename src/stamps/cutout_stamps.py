# Written by Anthony Pahl, 2/16/23

#This module is used to provide a way of using operating system dependent funcionality like reading or writing to a file system. Or just navigating a computer.
import os

#The subprocess module allows you to essentially run Unix commands from the script and store the output of those things in variables.
import subprocess

#Allows you to provide data to variables from the interpreter itself. I.e. with print(sys.argv[0]).

#You can then provide this data from the command-line like: python cutout_stamps.py "Hello World"
import sys

#This module is used to do high-level operations on files and the collection of files. I.e. shutil.copy('/home/user/Documents/config.ini', '/home/user/Backup/config.ini')


#shutil.copytree('source_dir', 'destination_dir') recursively copies an entire 

#You can use it to remove directories as well, shutil.rmtree('/home/user/Backup/old_backup')
import shutil

#Glob is used to find all the pathnames matching a specified pattern according to the rules used by the Unix shell. It's basically RegEx for files.

# Find all .txt files in the 'Documents' directory
#txt_files = glob.glob('/home/user/Documents/*.txt')
#print(txt_files)
import glob

#Used for creating and managing series and DataFrame files.

#Using dictionaries you can create 2D arrays with labels for columns and rows. I.e.:

#import pandas as pd

"""
data = 
{
    "Name": ["Alice", "Bob", "Charlie"],
    "Age": [25, 30, 35],
    "City": ["New York", "Los Angeles", "Chicago"]
}
df = pd.DataFrame(data)
print(df)
"""

#It can also read from a CSV file: df = pd.read_csv('path/to/file.csv')
#print(df.head())

#grouped = df.groupby('Name'.mean())
#print(grouped) <-- This would return a mean of all numerical values for a data set sorted by Name.

#ASCII and FITS are file formats, there isn't an ASCII or a FITS Python Package. Astropy, numpy, etc can all ready in
#these different files. Mainly just Astropy interacts with FITS files. They read in the contents and put it into their own
#proprietary system. An ASCII file would live on your computer, and then you'd use the other things in combination.

#Pandas is a general Python data analysis suite, astropy is more specific to astronomy. Astropy can read in FITS files.
#He uses the astropy input and if it's a table he'll turn it into a Pandas dataframe.

#They all have their own methods for creating 2D arrays. Pandas is the most sophisticated. Numpy is more barebones. Pandas
#also uses numpy under the good when it manipulates data. They're all necessary in different situations. 
import pandas as pd

#Allows you to do the pandas.series() and pandas.dataframe() functions without the pandas prefix.
from pandas import Series, DataFrame

#astropy is used to create and manipulate tables. I.e.
"""
data = {
    'name': ['Sirius', 'Betelgeuse', 'Rigel'],
    'distance': [8.6, 642.5, 863],
    'magnitude': [-1.46, 0.42, 0.12]
}
table = Table(data)
print(table)


  name    distance  magnitude
--------  ---------  ---------
  Sirius        8.6      -1.46
Betelgeuse     642.5      0.42
    Rigel       863       0.12
"""

#You can also define distances:
#from astropy import units as u
#distance = 10 * u.lightyear <--- 
#print(distance)


#output: 10.0 lightyear

#You can also convert from specified units 
#distance = 10 * u.lightyear

#distance_km = distance.to(u.kilometer)
#print(distance_km)

#output(distance_km)
#output: 9.460e+13 km

#Creating time objects:
#time = Time(''2023-07-01T00:00:00', format='isot', scale='utc')
#print(time)
#output: 2023-07-01 00:00:00.000

#astropy can also read FITS files:
#hdul = fits.open('example.fits')

#You can also access columns by their names and rows by their indicies

#print(table['name'])

#Access a row by its index
#print(table[0])

#?
#Why are we using astropy tables when we already have the tables from Pandas
from astropy.table import Table 

#Used to create arrays, you can also compute means with this, i.e. np.mean([1, 2, 3, 4, 5])

#It's also used for exponentiation, creating arrays, and trigonometric functions like np.sin(angle) or something else along those lines.
import numpy as np

#Copy can copy lists, i.e. shallow_copy = copy.copy(original_list)

#Changes made to the copied list will also reflect in the original list.
import copy

#astropy.io is a subpackage of astropy man to deal with input and output operations of data formats used in astronomy such as FITS or ASCII tables.

"""
from astropy.io import fits

with fits.open('example.fits' as hdul:')
    hdul.info()
    data = hdul[0].data
    header = hdul[0].header

#Print the data and header
print(data)
print(header)

from astropy.io import fits
import numpy as np

data = np.random.rand(100, 100)

# Create a HDU (Header/Data Unit)
hdu = fits.PrimaryHDU(data)

# Write the HDU to a new FITS file
hdu.writeto('new_example.fits', overwrite = True)

The header will be things like: 
[[0.1, 0.2, 0.3, 0.4, 0.5],
 [0.6, 0.7, 0.8, 0.9, 1.0],
 [1.1, 1.2, 1.3, 1.4, 1.5],
 [1.6, 1.7, 1.8, 1.9, 2.0],
 [2.1, 2.2, 2.3, 2.4, 2.5]]

 SIMPLE  =                    T / file does conform to FITS standard
 BITPIX  =                  -32 / number of bits per data pixel
 NAXIS   =                    2 / number of data axes

"""
from astropy.io import fits

#ASCII stands for American Standard Code for Information Interchange. In the context of data files, ASCII refers to plain text files where data is organized in a structured, readable format. An ASCII table is a text file where each row represents a record and each column represents a field.

#You can use ascii.read('example.csv') to read csv files.

#You can use ascii.write('csv') to write csv files.
#ascii.write(table, 'example_output.csv', overwrite=True)
from astropy.io import ascii

#The WCS (World Coordinate System) imports the WCS class from the astropy.WCS module. It is used to handle the mapping between pixel coordinates in an image and world coordinates, such as celestrial coordinates and other physical coordinates.

#Example of usage would be opening a FITS file and loading the WCS information. 

#Right Ascension (RA): Right ascension (RA) is the celestial equivalent of longitude.

#Declination (DEC): Declination (DEC) is the celestial sphere's equivalent of latitude.

#Definte world coordinates (RA, Dec)
#world_coords = SkyCoord([10.684, 10.685], [41.269, 41.270], unit=(u.deg, u.deg))

#Convert world coordinates into pixel coordinates:
#pixel_coords. = wcs.world_to_pixel(world_coords)

#It interfaces with the astropy UNIX module. WCS is something good to understand. It was good to have looked into this.
#FITS file:
#You have an array of pixel values
#You have a header with a bunch of different info: what instrument, how long the exposure was taken for, when it was taken, where in the sky the telescope was pointing when it took the image.
#It's important if we want to look at a specific galaxy you need to convert between the X, Y value in the image, and an actual coordinate in the sky. These SkyCoord objects can read in the WCS information in the header of a fits file.
#You can tell it I have this SkyCoord what coordinate does it correspond to. You can calculate it all by hand, but this makes it easier.
from astropy.wcs import WCS

#Imports the SkyCord class from the astropy.coordinates submodule.

#coord = SkyCoord(ra=10.684*u.deg, dec=41.269*u.deg, frame='icrs')

#This creates a SkyCoord object with right ascension (RA) of 10.684 degrees and declination (Dec) of 41.269, for example.

#You can then calculate the angular seperation of two instances of SkyCoord with things like coord1.separation(coord2)
#WCS is the coordinate system of the thing you're looking at while SkyCoord is one specific coordinate. SkyCoord is just like one coordinate.
from astropy.coordinates import SkyCoord

#Imports the Time class from the astropy.time submodule. You can to arthmetic to these objects or even use things like terrestial time.

#Times for different things in astronomy. There are different systems that people use and some of them are fairly obtubse. The nice thing about high redshift galaxies is that they don't change very fast. The image of a super distant galaxy and 20 years in the future is going to look almost exactly the same. It's not that important on Earth when these things were taken.
from astropy.time import Time

#Imports the Cutout2D class from the astropy.nddata module. This is used for the actual cutting out of images from larger 2D arrays, such as images.

#It maintains World Coordinate System (WCS) information for the cutout, ensuring spatial information is preserved.

#Basic usage looks something like this:
#"cutout = Cutout2D(data, position, size, wcs=wcs)"

#This is the backbone of the code. It's basically doing the stamp cutting for us. Everything else ensures this works correctly. It takes the data as an input. It also could take like a SkyCoord coordinate, decides the stamp of what you want, and then takes it.
from astropy.nddata import Cutout2D

#Imports astropy.units as u allowing you to do multiplication to numbers to assign units.
#If you multiply a regular floating point number by an astropy unit it's now a quantity object. Once you multiply
#it by a unit it's no longer a normal number. 
import astropy.units as u

#Imports the  skycoord_to_pixel function from the astropy.wcs.utils submodule, allowing you to perform operations to SkyCoord objects. 

#Converts the celestial coordinates (e.g., RA and Dec) to pixel coordinates in an image.

#After creating a SkyCoord object from say opening a Flexible Image Transport System (FITS) file, it allows you to convert Right Ascension, Declination coordinates to pixel coordinates

#Example: pixel_coords = skycoord_to_pixel(sky_coord, wcs)
from astropy.wcs.utils import skycoord_to_pixel 

#This is a Python wrapper for the SExtractor (Source Extractor) software. Source Extractor is a popular program used in astronomy for detecting and measuring objects in astronomical images.


#Basically you specify the things you want to extract rom the FITS file:
#sew = sewpy.SEW(params=["X_IMAGE", "Y_IMAGE", "MAG_AUTO"], config={"DETECT_MINAREA": 5, "DETECT_THRESH": 1.5})

#Then you can get the results of those things with result = sew(fits_image_path)

#Notes in Google Docs
#import sewpy

#Imports all variables from vars
from vars import *


#Defines a function named "invert_sqrt_fits" that takes to arguments. The input FITS file name (path), and the arr_out "the output FITS file name path"
def invert_sqrt_fits(arr_in,arr_out):
    #Open the fits file from the input path.
    with fits.open(arr_in) as hdu:
        #Data is set to that of the primary HDU (Header + Data Units) of the FITS file.
        #We're reading in the scientific image. So later in the code he does this on the weight image. The primary HDU is the weighted information.

        #You're downloading the science mosaic and the weight mosaic. The weight mosaic is the same dimensions as the science but at every pixel value it's telling you the error for the science image.

        #The weight if converted to the sigma then then the sigma is the amount of standard deviations or percision of the pixel value. GalFit just wants this as a sigma image.

        #Pixels with higher weight are more confident. There's less of an intreval of what it could be.

        #It's just the WEIGHT image not the WEIGHTED image. 
        data = hdu[0].data
        #Converts the weight image to the sigma image.
        data2 = np.sqrt(1 / data)
        #We then set the primary HDU equal to that of data 2.
        hdu[0].data = data2

        #After that we write the data object we opened to the output file. We're basically reformatting the data by putting it all over the square root of one.

        #He reads in the file into memory then he's writing it out to a new file name.
        hdu.writeto(arr_out,overwrite=True)

#We are reading pickled Pandas data from '{}/cat/ceers_radec.pkl'.format(res_path)' which uses Python string formatting to create a file path.

#res_path is from vars and contains the base directory where the file in question is located.

#This pickled Pandas data frame is then assigned to pcrs
pcrs = pd.read_pickle('/Users/rschmidt-eder/cassi2024/data/jc_xi_ion_prop.pkl')

pcrs = pcrs.loc[pcrs.index.str.contains('JADES')]

# This code snippet is used to process and store image data, headers, and World Coordinate System (WCS) information from a set of FITS files into seperate dictionaries.

#Deep copies (files including metadata are copied over from the science image directory [I think]) for the image, header, and world coordinate system directory.

####
#It copies the dictionaries, and then instead of the value being the file name it’s the mosaic read into memory.
#^ This entire thing returns a mosaic read into memory. Otherwise you'd have to read in the mosaic every time. It'll slow down your code by a lot if you don't read them in first. 
####

sci_img_dir = copy.deepcopy(sci_dir)
sci_header_dir = copy.deepcopy(sci_dir)
sci_wcs_dir = copy.deepcopy(sci_dir)

#The outer loop iterates over the keys of sci_dir. Each fi represents a filter in the directory structure.
for fi in sci_dir:
    #Iterates over the fields within each filter. And for each field (I'm guessing this means image) it opens the filter and field as variable hdu.
    for field in sci_dir[fi]:

        #The header, the WVS, and the img all go into their own dictionary to have for later #####
        with fits.open(sci_dir[fi][field]) as hdu:
            #The image data from the first HDU is extracted.
            img = hdu[0].data
            #A WCS object from the header of the first HDU.
            wcs = WCS(hdu[0].header)
            #Extracts the header from the first HDU.
            header = hdu[0].header
        #Stores image data in the science image directory under the same filter and field keys.
        sci_img_dir[fi][field] = img
        #Stores the header in sci_header_dir under the same filter and field keys.
        sci_header_dir[fi][field] = header
        #Stores the WCS object in 'sci_wcs_dir' under the same filter and field keys.
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
#Commented out for now. This is for segment images.
"""
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
"""
# Determine image membership

#Figures out the field each object lives in. It makes a sky coordinate for that object and sees if that sky coordinate is in the image. If it is then it puts that into a new column in the Pandas dataframe and it calls up field_filter. The value would be an actual field where the object lays in.
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
pcrs.loc[:, ['field_F115W']].to_pickle('{}/stamps/field_mem.pkl'.format(pdata_path))

# Cutout stamps
#It's a for loop but the for loop is itereating over every galaxy in the catalog.
#Depending on the filter he needs the mosaic to do cutouts of it.
#He could have it in the for loop but then for every 100+ galaxies he'd have to read in that GB file every time.
#^^ Doing it intially prevents you from having to do it every time in the for loop. He puts it into a dictionary it's easier to call it later.
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
            

