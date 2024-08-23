import matplotlib.pyplot as plt
import pandas as pd

# Effective radii in kpc and xi_ion values
effective_radii_kpc = [0.104, 0.374, 1.410, 0.217, 0.293, 0.225, 0.109, 0.413, 0.271, 0.632, 0.171, 2.410, 2.272, 0.477]

folders = [
    "JADES-4297", "JADES-4404", "JADES-5329", "JADES-6246",
    "JADES-9452", "JADES-16625", "JADES-17260", "JADES-17777",
    "JADES-18846", "JADES-18970", "JADES-19342", "JADES-19519", "JADES-10000626",
    "JADES-10013618"
]

pjades = pd.read_pickle("/Users/rschmidt-eder/cassi2024/data/jc_xi_ion_prop.pkl")

xi_ions = pjades.loc[folders, 'xiion']

# Create scatter plot
plt.figure(figsize=(10, 6))
plt.scatter(effective_radii_kpc, xi_ions, color='blue')
plt.title('Scatter Plot of xi_ion vs. Effective Radii (kpc)', fontsize=20)
plt.xlabel('Effective Radii (kpc)', fontsize=16)
plt.ylabel('log10(xi_ion / Hz erg^-1)', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)
plt.show()
