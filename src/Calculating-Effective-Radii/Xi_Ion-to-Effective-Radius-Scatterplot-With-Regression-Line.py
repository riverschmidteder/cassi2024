import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
import pandas as pd
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import os

# Cosmology model
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

simplemorphologygalaxies = [
    "JADES-4297", "JADES-4404", "JADES-5329", "JADES-6246",
    "JADES-9452", "JADES-16625", "JADES-17260", "JADES-17777",
    "JADES-18846", "JADES-18970", "JADES-19342", "JADES-19519", "JADES-10000626",
    "JADES-10013618"
]

base_dir = "/Users/rschmidt-eder/cassi2024/processed_data/galfit"

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


#Have the iPython commands directly in here.

#Kpc values will need to be regenerated from the 

#Should've connected to the pickle file from dconv.

# Provided kpc values and xi_ion values
kpc_values = np.genfromtxt("Sizes-in-Kpcs.txt")
psizes = pd.read_pickle("/Users/rschmidt-eder/cassi2024/results/size/isophot_size_deconv.pkl")
pproperties = pd.read_pickle("/Users/rschmidt-eder/cassi2024/data/jc_xi_ion_prop.pkl")
xi_ions = pproperties.loc[psizes.index, "xiion"].values
xi_ionlow = pproperties.loc[psizes.index, "xiion_lo"].values
xi_ionhigh = pproperties.loc[psizes.index, "xiion_hi"].values
galaxynames = psizes.index.values
redshifts = pproperties.loc[psizes.index, "z"].values

kpc_values_error = np.empty(len(kpc_values))
kpc_values_error[:] = np.nan

for i, folder in enumerate(galaxynames):
    if folder in simplemorphologygalaxies:
        # Construct the path to the FITS file
        fits_file = os.path.join(base_dir, folder, "F115W_single", f"{folder}.F115W_single.fit.fits")
        
        # Open the FITS file
        with fits.open(fits_file, ignore_missing_simple=True) as hdul:
            header = hdul[2].header
            image_data = hdul[1].data

        # Extract parameters and uncertainties from the header
        params = extract_params_from_header(header)

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
        sRe_kpc = cosmo.kpc_proper_per_arcmin(redshifts[i]).value * sRe_corr_arcsec * u.arcsec.to(u.arcmin)
        Re_kpc = cosmo.kpc_proper_per_arcmin(redshifts[i]).value * Re_corr_arcsec * u.arcsec.to(u.arcmin)

        kpc_values_error[i] = sRe_kpc

# Remove NaN values for both arrays
mask = ~np.isnan(kpc_values) & ~np.isnan(xi_ions)
kpc_values_clean = np.array(kpc_values)[mask]
xi_ions_clean = np.array(xi_ions)[mask]

# Calculate Spearman correlation
sr, sp = spearmanr(kpc_values_clean, xi_ions_clean)
print("SPEARMAN: " + str(sr), str(sp))

# Fit linear regression via least squares with numpy.polyfit
b, a = np.polyfit(kpc_values_clean, xi_ions_clean, deg=1)

# Print slope and intercept
print(f"Slope (b): {b}")
print(f"Intercept (a): {a}")

# Initialize layout
fig, ax = plt.subplots(figsize=(10, 6))

# Add scatterplot
ax.errorbar(kpc_values, xi_ions, xerr=kpc_values_error, yerr=[xi_ionlow, xi_ionhigh], ms=2, alpha=0.7, color="k", linestyle="", marker="o")

# Create sequence of 100 numbers from min(x) to max(x)
xseq = np.linspace(min(kpc_values_clean), max(kpc_values_clean), num=100)

# Plot regression line
ax.plot(xseq, a + b * xseq, color="k", lw=2.5)

# Add titles and labels with larger font sizes
ax.set_title('Scatter Plot of xi_ion vs. Effective Radii (kpc) with Regression Line', fontsize=20)
ax.set_xlabel('Effective Radii (kpc)', fontsize=16)
ax.set_ylabel('Ionizing Photon Production Efficiencies', fontsize=16)
ax.tick_params(axis='both', which='major', labelsize=14)

# Show grid
ax.grid(True)

plt.show()
