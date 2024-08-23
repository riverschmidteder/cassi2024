import matplotlib.pyplot as plt

# Effective radii measurements in kpc
effective_radii_kpc = [0.104, 0.374, 1.410, 0.217, 0.293, 0.225, 0.109, 0.413, 0.271, 0.632, 0.171, 2.410, 2.272, 0.477]

# Create histogram
plt.figure(figsize=(10, 6))
plt.hist(effective_radii_kpc, bins=10, edgecolor='black', color='orange')
plt.title('Histogram of Effective Radii in kpc')
plt.xlabel('Effective Radii (kpc)')
plt.ylabel('Frequency')
plt.grid(True)
plt.show()