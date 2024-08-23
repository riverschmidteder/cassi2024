import os
import subprocess

# Base directory where all your galaxy folders are located
base_dir = "/Users/rschmidt-eder/cassi2024/processed_data/stamps"

# Folders to exclude
exclude_folders = [
    "JADES-3334", "JADES-6002", "JADES-22251", "JADES-22924", "JADES-10011849"
]

# Loop through all items in the base directory
for folder in os.listdir(base_dir):
    # Skip excluded folders
    if folder in exclude_folders:
        continue

    folder_path = os.path.join(base_dir, folder, "F115W")
    
    # Ensure it's a directory
    if os.path.isdir(folder_path):
        # Look for the FITS file in the directory
        for file in os.listdir(folder_path):
            if file.endswith("_F115W_fit.fits"):
                fits_file_path = os.path.join(folder_path, file)
                print(f"Attempting to open: {fits_file_path}")
                # Directly use the full ds9 command
                try:
                    result = subprocess.run(
                        f"ds9 -lock frame physical -lock scalelimits yes -lock colorbar yes -cmap sls {fits_file_path}[1] {fits_file_path}[2] {fits_file_path}[3] -frame first -match scalelimits -zscale ./stamps/*bpm.2as.fits &",
                        shell=True,
                        check=True,
                        capture_output=True
                    )
                    print(f"Command executed successfully: {result}")
                except subprocess.CalledProcessError as e:
                    print(f"Error occurred: {e}")
