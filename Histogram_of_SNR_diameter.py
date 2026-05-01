# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 12:21:58 2026

@author: lily zaidi
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 
# read histogram text file
df = pd.read_csv("histogram_snr_data.txt", sep="|", names=["ID","Ratio","Diam_pc"])

# remove spaces and convert to numbers
df = df.apply(lambda x: x.astype(str).str.strip()).astype(float)

# keep only real SNRs (Ratio > 0.4)
snr_df = df[df["Ratio"] >= 0.4]

# plot histogram of diameter
plt.hist(snr_df["Diam_pc"], bins=10, edgecolor="red", color='orange')

plt.xlabel("Diameter (pc)")
plt.ylabel("Number of SNRs")
plt.grid()
plt.title("SNR Diameter Distribution (Ratio > 0.4)")


plt.show()

hii_df = df[df["Ratio"] < 0.4]

plt.hist(hii_df["Diam_pc"], bins=10, edgecolor="purple", color='magenta')
plt.xlabel("Diameter (pc)")
plt.grid()
plt.ylabel("Number of H II Regions")
plt.title("H II Region Diameter Distribution")

plt.show()


plt.hist(snr_df["Diam_pc"], bins=10, alpha=0.6, color='orange',label="SNR (Ratio > 0.4)")
plt.hist(hii_df["Diam_pc"], bins=10, alpha=0.6, color='magenta',label="HII (Ratio ≤ 0.4)")

plt.xlabel("Diameter (pc)")
plt.ylabel("Number of Objects")
plt.grid()
plt.title("Diameter Distribution: SNR vs HII")
plt.legend()
plt.show()
plt

print("Number of HII regions: ", len(hii_df))
print("Number of SNRs: ", len(snr_df))


#representing the ratio as a cut off point, 
#anything under 0.4 is a HII region photoionized 
#anything <= to 0.4 is a SNR region 
import pandas as pd
import matplotlib.pyplot as plt

snr_df = df[df["Ratio"] > 0.4]
hii_df = df[df["Ratio"] < 0.4]

# plot
plt.figure(figsize=(8,6))

plt.scatter(snr_df["Diam_pc"], snr_df["Ratio"],
            color="red", label="SNR (Ratio ≥ 0.4)", s=25)

plt.scatter(hii_df["Diam_pc"], hii_df["Ratio"],
            color="blue", label="HII (Ratio < 0.4)", s=25, alpha=0.7)

# threshold line
plt.axhline(0.4, linestyle="--", color="black", linewidth=2,
            label="Matthewson-Clarke threshold (0.4)")

# labels and title
plt.xlabel("Diameter (pc)")
plt.ylabel("Matthewson-Clarke Ratio")
plt.title("Diameter vs Ratio for SNR and HII Regions")

plt.legend()
plt.grid(alpha=0.3)


#rep_diam = 39.5249
#x_err_value = rep_diam * 0.287

"""x_pos = 39.5249
y_pos = 0.309

plt.errorbar(x_pos, y_pos, xerr=x_err_value, fmt='o', color='black', 
             capsize=5, elinewidth=2, label="Mean Spatial Uncertainty")
plt.text(x_pos, y_pos + 0.15, '±28.7% Diameter Error', 
         ha='center', fontsize=9, fontweight='bold')

plt.show()"""

#culmative dist plot 
import numpy as np

plt.figure(figsize=(8,6))

#frequency of ratio 
plt.figure(figsize=(8,6))
plt.hist(df["Ratio"], bins=20, alpha=0.5, edgecolor="seagreen", color="lime")
plt.xlabel("Ratio")
plt.ylabel("Frequency")
plt.title("Frequency Distribution of Matthewson-Clarke Ratio")
plt.axvline(0.4, linestyle="--", color="blue", linewidth=2, label="Threshold (0.4)")
plt.legend()
plt.grid(alpha=0.3)
plt.show()

median_sii = np.median(snr_df)
print(median_sii)


