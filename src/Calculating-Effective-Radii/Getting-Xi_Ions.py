import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import pandas as pd
import pickle

# Directory containing the folders with FITS files
base_dir = "/Users/rschmidt-eder/cassi2024/processed_data/galfit"

# Path to the pickle file with xiion values
pickle_file_path = "/Users/rschmidt-eder/cassi2024/data/jc_xi_ion_prop.pkl"

# List of folders
folders = [
    "JADES-4297", "JADES-4404", "JADES-5329", "JADES-6246",
    "JADES-9452", "JADES-16625", "JADES-17260", "JADES-17777",
    "JADES-18846", "JADES-18970", "JADES-19342", "JADES-19519", "JADES-10000626",
    "JADES-10013618"
]

effective_radii_kpc = [0.104, 0.374, 1.410, 0.217, 0.293, 0.225, 0.109, 0.413, 0.271, 0.632, 0.171, 2.410, 2.272, 0.477]

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

# Load xiion values from the pickle file
with open(pickle_file_path, 'rb') as file:
    xi_ion_data = pickle.load(file)

# Cosmology model
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

# Create a figure with subplots, 2 rows 7 columns
fig, axes = plt.subplots(2, 7, figsize=(16, 8))
axes = axes.flatten()

xi_ion_values = []

for i, folder in enumerate(folders):
    try:
        # Construct the path to the FITS file
        fits_file = os.path.join(base_dir, folder, "F115W_single", f"{folder}.F115W_single.fit.fits")
        
        # Open the FITS file
        with fits.open(fits_file, ignore_missing_simple=True) as hdul:
            header = hdul[2].header
            image_data = hdul[1].data

        # Extract parameters and uncertainties from the header
        params = extract_params_from_header(header)

        # If the galaxy name is in the pickle file with the xiion
        if folder in xi_ion_data:
            xi_ion = xi_ion_data[folder]['xiion']
            xi_ion_values.append(xi_ion)
        else:
            raise ValueError(f"xiion value for {folder} not found in pickle file")

        # Calculate the effective radius in arcseconds
        Re = params['Re']
        AR = params['AR']
        sRe = params['sRe']
        sAR = params['sAR']

        Re_corr = Re * np.sqrt(AR)
        sRe_corr = Re_corr * np.sqrt((sRe / Re)**2 + (1/2)*(sAR / AR)**2)
        Re_corr_arcsec = Re_corr * 0.03
        sRe_corr_arcsec = sRe_corr * 0.03
        size_arcsec = np.pi * Re_corr_arcsec**2
        Re_kpc = cosmo.kpc_proper_per_arcmin(xi_ion).value * Re_corr_arcsec * u.arcsec.to(u.arcmin)

        # Plot the image
        ax = axes[i]
        ax.imshow(image_data, cmap='gray', origin='lower')
        title = f"{folder}\nRe: {Re_corr_arcsec:.3f} (arcsec), \n{Re_kpc:.3f} (kpc)"
        ax.set_title(title)
        ax.axis('off')

    except Exception as e:
        print(f"Error processing file {fits_file}: {e}")

# Adjust layout and show the plot
plt.tight_layout()
plt.show()

# Pair the kpc values with the galaxies and xi_ion values
paired_values = list(zip(folders, effective_radii_kpc, xi_ion_values))

# Display the paired values
for folder, kpc, xi_ion in paired_values:
    print(f"{folder}: kpc={kpc}, xi_ion={xi_ion}")
