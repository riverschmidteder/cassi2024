from astropy.io import fits

hdul = fits.open("/Users/rschmidt-eder/cassi2024/processed_data/galfit/JADES-4297/F115W_single/JADES-4297.F115W_single.fit.fits")

hdr = hdul[2].header

print(hdr)
