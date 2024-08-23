import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u

# Initialize the cosmology model
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

psizes = pd.read_pickle("/Users/rschmidt-eder/cassi2024/results/size/isophot_size_deconv.pkl")
pproperties = pd.read_pickle("/Users/rschmidt-eder/cassi2024/data/jc_xi_ion_prop.pkl")
# Provided Re_arcsecond values and their corresponding redshifts
"""
Re_arcsecs = [
    0.25272746577723526, 0.10854520185920347, 0.04980395428531486,
    0.22549699372785773, np.nan, 0.08804178451262737,
    0.017608356902525477, 0.0704334276101019, 0.48061461009940576,
    0.21203269501889063, 0.04658733336055001, 0.15945045340214858,
    0.030498568792980265, 0.5332164021143376, 0.06588443867328636,
    0.04658733336055001, 0.08626297924000305, 0.05282507070757643,
    0.09317466672110002, 0.3140019841219731, 0.0393734830061246,
    0.15847521212272928, 0.19765331601985908, 0.09317466672110002,
    0.05840031302034161, 0.1545127047401513, np.nan, np.nan,
    0.09644494276114998, 0.16138325673924228, np.nan, 0.15847521212272928,
    0.06588443867328636, 0.030498568792980265, 0.06819687303755759,
    0.19607836409974674, np.nan, 0.043131489620001524,
    0.05568251366512799, 0.017608356902525477, 0.12450988571328714,
    0.2163751121405006, 0.0393734830061246, 0.07260111540268635,
    0.05282507070757643
]
"""

Re_arcsecs = psizes.RE_F115W_single.values
redshifts = pproperties.loc[psizes.index, "z"].values


"""
redshifts = [
    4.46240056, 2.98893705, 3.5975703, 2.57535419, 2.69237274,
    5.94091625, 5.9143494, 5.07373332, 3.92578392, 3.58744148,
    6.62780035, 5.56371609, 4.88346005, 6.0987044, 4.1314925,
    4.77275453, 3.1497631, 6.33219992, 3.72292189, 6.32447915,
    5.97088292, 3.31859997, 3.60271249, 5.88603402, 2.84292996,
    3.0860869, 5.79547999, 3.70106334, 3.0131326, 3.46567939,
    6.70318237, 2.8054645, 4.02046921, 6.70957761, 5.76013789,
    3.58505689, 5.93408965, 5.55897135, 5.61276155, 4.04177546,
    4.14634457, 3.57610463, 4.22654442, 4.8037697, 5.12061427
]
"""

# Create a list to store the calculated Re_kpc values
Re_kpc_values = []

# Loop through each galaxy's Re_arcsecond and redshift value
for Re_arcsec, z in zip(Re_arcsecs, redshifts):
    if not np.isnan(Re_arcsec) and not np.isnan(z):
        Re_kpc = cosmo.kpc_proper_per_arcmin(z).value * Re_arcsec * u.arcsec.to(u.arcmin)
        Re_kpc_values.append(Re_kpc)
    else:
        Re_kpc_values.append(np.nan)

# Display the results
for i, (Re_arcsec, z, Re_kpc) in enumerate(zip(Re_arcsecs, redshifts, Re_kpc_values)):
    print(f"Galaxy {i+1}: Re_arcsec = {Re_arcsec:.6f}, Redshift = {z:.6f}, Re_kpc = {Re_kpc:.6f} kpc")

# Optionally plot the results
plt.figure(figsize=(10, 6))
plt.plot(redshifts, Re_kpc_values, 'o', label='Galaxy Sizes in kpc')
plt.xlabel('Redshift')
plt.ylabel('Effective Radius (kpc)')
plt.title('Galaxy Sizes Converted to Kiloparsecs')
plt.grid(True)
plt.legend()
plt.show()
# Print the array of Re_kpc values
print("Array of Re_kpc values (in kpc):", np.array(Re_kpc_values))

with open("Sizes-in-Kpcs.txt", "w") as outfile:
    outfile.write("\n".join(str(item) for item in Re_kpc_values))

