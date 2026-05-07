import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from astropy.io import fits
from astropy.units import Unit

print("hello world") #to test if the code is running correctly


def load_fits(fname):
    with fits.open(fname) as hdu:
        data =  hdu[1].data  
    return pd.DataFrame(data)

df= load_fits('/Users/emmalinares/Downloads/hlsp_jades_jwst_nirspec_goods-s_prism-line-fluxes_v1.1_catalog.fits') #load the data using the function we just created
print(df.columns) #print the names of the columns in the data to understand what information we have
print(df.shape) #print the shape of the data to see how many galaxies we have and how many columns of information we have for each galaxy

# we want 3 things: Galaxy identification , Redshift, and spectra!

columns = [
    'NIRSpec_ID', #galaxy id
    'RA_TARG', #right ascension of galaxy
    'Dec_TARG', #declination of galaxy 
    'z_Spec', #redshift of galaxy
    'z_Spec_flag', #redshift quality flag
    'HA_6563_flux', #flux of the H-alpha line -star formation indicator 
    'HB_4861_flux', #flux of the H-beta line  -star formation indicator
    'O3_5007_flux', #flux of the OIII 5007 line -AGN indicator
    'O3_4959_flux', #flux of the OIII 4959 line -AGN indicator
]

df= df[columns] #keep only the columns we are interested in
df= df.dropna()

df = df[df["z_Spec_flag"].isin(["A", "B", "C"])]
# ^^ keepings only the secure redshift galaxies w high quality flags
df = df[df["z_Spec"] > 5] 

print('shape')

plt.plot(df['z_Spec'], df['HA_6563_flux'], 'o') #plot the redshift vs the H-alpha flux to see if there is any relationship between them
plt.xlabel('Redshift')
plt.ylabel('H-alpha Flux')
plt.title('star Formation in Early Galaxies')
plt.savefig('figures/ha_flux_vs_redshift.png', dpi = 300, bbox_inches = "tight")
plt.show()

df= df[
    (df["HB_4861_flux"] > 0) &
    (df["O3_5007_flux"] > 0)
]

df["O3_Hb_ratio"] = (
    df["O3_5007_flux"] /
    df["HB_4861_flux"]
)

plt.scatter(df["z_Spec"], df["O3_Hb_ratio"], color = "red") #plot the redshift vs the OIII/H-beta ratio to see if there is any relationship between them
plt.yscale("log")
 #set the y-axis to a logarithmic scale bc emission line ratios
# span several orders of magnitude
plt.xlabel('Redshift')
plt.ylabel('[OIII] /H-beta Ratio')
plt.title('AGN Diagnostics in Early Galaxies')
plt.savefig("figures/o3_hb_ratio_vs_redshift.png", dpi=300, bbox_inches = "tight")
plt.show()

print(df.shape) #print the shape of the data after filtering to see how many galaxies we have left

#finding AGN candidates using the OIII/H-beta ratio. 

outlier = df.sort_values(
    by= "O3_Hb_ratio",
    ascending = False
)

print(outlier.head()) 

highest_z = df.sort_values(
    by= "z_Spec",
    ascending = False
)

print(highest_z.head())

#HOORAY! we have our candidates for AGN 

candidate_1= df[df["NIRSpec_ID"] == 20049790] 
print(candidate_1.T)

# so far appears to be a possible extreme AGN , with extremely
#high OIII/H-beta ratio and a redshift of 6.5, but could be due to noisy data, 
# and have weak H-beta so.. ??? 

candidate_2 = df[df["NIRSpec_ID"] == 13609 ]
print(candidate_2.T) 

#high redshift of z= 6.94
#stong Ha and Hb flux  emission, indicates active star formation
# 03/hb ratio of 7.8 which is very high, suggests strong ionization
#  good AGN candidate! 


