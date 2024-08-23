import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import pandas as pd

# Directory containing the folders with FITS files
base_dir = "/Users/rschmidt-eder/cassi2024/processed_data/galfit"

# Path to the pickle file with z_spec values
pickle_file_path = "/Users/rschmidt-eder/cassi2024/data/jc_xi_ion_prop.pkl"

# List of folders
folders = [
    "JADES-4297", "JADES-4404", "JADES-5329", "JADES-6246",
    "JADES-9452", "JADES-16625", "JADES-17260", "JADES-17777",
    "JADES-18846", "JADES-18970", "JADES-19342", "JADES-19519", "JADES-10000626",
    "JADES-10013618"
]

# Function to extract parameters and their uncertainties from the FITS file header
def extract_params_from_header(header):
    params = {}
    try:
        re_value, re_uncertainty = header['1_RE'].split(' +/- ')
        ar_value, ar_uncertainty = header['1_AR'].split(' +/- ')
        
        params['Re'] = float(re_value)
        params['sRe'] = float(re_uncertainty)
        params['AR'] = float(ar_value)
        params['sAR'] = float(ar_uncertainty)
    except KeyError as e:
        raise ValueError(f"Required parameter {e} not found in header")
    return params

# Load z_spec values from the pickle file
z_spec_data = pd.read_pickle(pickle_file_path)

# Cosmology model
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

# Create a list to store the calculated Re_kpc and corresponding folder names
galaxies_data = []

for folder in folders:
    try:
        fits_file = os.path.join(base_dir, folder, "F115W_single", f"{folder}.F115W_single.fit.fits")
        print(f"Processing file: {fits_file}")
        
        with fits.open(fits_file, ignore_missing_simple=True) as hdul:
            header = hdul[2].header
            image_data = hdul[1].data

        params = extract_params_from_header(header)
        
        if folder in z_spec_data.index:
            z_spec = z_spec_data.loc[folder, 'z']
        else:
            raise ValueError(f"z_spec value for {folder} not found in pickle file")

        Re = params['Re']
        AR = params['AR']
        sRe = params['sRe']
        sAR = params['sAR']

        Re_corr = Re * np.sqrt(AR)
        sRe_corr = Re_corr * np.sqrt((sRe / Re)**2 + (1/2)*(sAR / AR)**2)
        Re_corr_arcsec = Re_corr * 0.03
        sRe_corr_arcsec = sRe_corr * 0.03
        size_arcsec = np.pi * Re_corr_arcsec**2
        Re_kpc = cosmo.kpc_proper_per_arcmin(z_spec).value * Re_corr_arcsec * u.arcsec.to(u.arcmin)
        
        galaxies_data.append((Re_kpc, image_data, folder, Re_corr_arcsec))

    except Exception as e:
        print(f"Error processing file {fits_file}: {e}")

# Sort the galaxies data by Re_kpc in ascending order
galaxies_data.sort(key=lambda x: x[0])

# Create a figure with subplots
fig, axes = plt.subplots(2, 7, figsize=(16, 8))
axes = axes.flatten()

for i, (Re_kpc, image_data, folder, Re_corr_arcsec) in enumerate(galaxies_data):
    ax = axes[i]
    ax.imshow(image_data, cmap='gray', origin='lower')
    title = f"{folder}\nRe: {Re_corr_arcsec:.3f} (arcsec), \n{Re_kpc:.3f} (kpc)"
    ax.set_title(title)
    ax.axis('off')

plt.tight_layout()
plt.show()


