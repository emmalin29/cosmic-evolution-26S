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

# we want 3 things: Galaxy identification , Redshift, and spectra!

columns = [
    'NIRSpec_ID', #galaxy id
    'RA_TARG', #right ascension of galaxy
    'Dec_TARG', #decination of galaxy 
    'z_Spec', #redshift of galaxy
    'z_Spec_flag', #redshift quality flag
    'HA_6563_flux', #flux of the H-alpha line -star formation indicator 
    'HB_4861_flux', #flux of the H-beta line  -star formation indicator
    'O3_5007d_flux', #flux of the OIII 5007 line -AGN indicator
    'O3_4959_flux', #flux of the OIII 4959 line -AGN indicator
]

print('shape')
print(df.shape)
print(df.columns)

plt.plot(df['z_Spec'], df['HA_6563_flux'], 'o') #plot the redshift vs the H-alpha flux to see if there is any relationship between them
plt.xlabel('Redshift')
plt.ylabel('H-alpha Flux')
plt.title('Galaxy Redshift vs H-alpha Flux')
plt.show()

