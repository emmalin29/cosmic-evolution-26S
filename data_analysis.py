import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from astropy.io import fits
from astropy.units import Unit
from astroquery.mast import Observations 
from astropy.table import Table
import os

pd.set_option('display.max_columns', None) #set pandas option to display all columns when we print the data


def load_fits(fname):
    with fits.open(fname) as hdu:
        data =  hdu[1].data  
        header= hdu[0].header
        print(header)
    return pd.DataFrame(data)

df= load_fits('/Users/emmalinares/Downloads/hlsp_jades_jwst_nirspec_goods-s_prism-line-fluxes_v1.1_catalog.fits') #load the data using the function we just created

print(df.columns) #print the names of the columns in the data to understand what information we have
print(df.shape) #print the shape of the data to see how many galaxies we have and how many columns of information we have for each galaxy

# we want 3 things: Galaxy identification , Redshift, and spectra!

columns = [
    'NIRSpec_ID', #galaxy id
    'NIRCam_ID', #id for the viewer images
    'RA_TARG', #right ascension of galaxy
    'Dec_TARG', #declination of galaxy 
    'z_Spec', #redshift of galaxy
    'z_Spec_flag', #redshift quality flag.   # star forming indicators below: 
    'HA_6563_flux', #H-alpha line 
    'HA_6563_err', #error on the H-alpha 
    'HB_4861_flux', #H-beta line 
    'HB_4861_err', #error on the H-beta 
    'O2_3727_flux', #OII 3727 line
    'O2_3727_err', #error on the OII 3727      #agn indicators below:
    'O3_5007_flux', # OIII 5007 line 
    'O3_5007_err', #error on the OIII 5007 
    'Blnd_He2_O3_1650_flux', #flux of the blended HeII 1640 line
    'Blnd_He2_O3_1650_err', #error on the blended HeII 1640   #uhmmm perchance needed below:
    'O3_4959_flux', #flux of the OIII 4959 line -AGN indicator
    'O3_4959_err', #error on the OIII 4959 flux measurement

]
#future use for JADES viewer pics (potentially)
print(df['NIRCam_ID'][0:10])

df= df[columns] #keep only the columns we are want
df = df[df["z_Spec_flag"].isin(["A", "B", "C"])]
# ^^ keeping only the secure redshift galaxies w high quality flags
df = df[df["z_Spec"] > 5] 

#print('shape')

sf_df= df[
    (df["HA_6563_flux"] > 0)
]
df["O3_Hb_ratio"] = df["O3_5007_flux"] / df["HB_4861_flux"]

sfsf_df= df[
    (df["HA_6563_flux"] >0) &
    (df["O3_Hb_ratio"] < 5)
]

agnagn_df = df[
    (df["O3_Hb_ratio"] > 5)
]

agn_df= df[
    (df["HB_4861_flux"] > 0) &
    (df["O3_5007_flux"] > 0)
]

agn_df["O3_Hb_ratio"] = (
    df["O3_5007_flux"] /
    df["HB_4861_flux"]
)

#NEW IDEA GONNA INCLUDE REDSHIFT ERROR muahahah
flag_markers = {
    "A": "o",
    "B": "s",
    "C": "^"
}
plt.figure(figsize=(8,6))

for flag, marker in flag_markers.items():
    temp = df[df["z_Spec_flag"] == flag]


    plt.scatter(temp["z_Spec"],temp["O3_Hb_ratio"],color="gray",marker=marker, alpha=0.4,s=35,label=f"Flag {flag}")
plt.scatter(sfsf_df["z_Spec"], sfsf_df['HA_6563_flux'],marker = "*", color = "blue", label = "Star forming galaxy candidates")
plt.scatter(agnagn_df["z_Spec"], agnagn_df["O3_Hb_ratio"], marker = "*", color = "red", label = "AGN galaxy candidates") #plot the redshift vs the OIII/H-beta ratio to see if there is any relationship between them
plt.yscale("log")
plt.xlim(5, 9.5)
plt.xlabel('Redshift')
plt.ylabel(r"$\frac{[O_{III}]}{H_{\beta}}$ Ratio")
plt.title("Star Formation and AGN Diagnostics in Early Galaxies (High Redshift )")
plt.legend()
plt.show()


# these plots below are just the sepearte star forming and agn plots idk if i still need

#star forming galaxy plot 
#plt.plot(sf_df['z_Spec'], sf_df['HA_6563_flux'], 'o') #plot the redshift vs the H-alpha flux to see if there is any relationship between them
#plt.xlabel('Redshift')
#plt.ylabel('H-alpha Flux')
#plt.title('star Formation in Early Galaxies')
#plt.savefig('figures/ha_flux_vs_redshift.png', dpi = 300, bbox_inches = "tight")
#plt.show()

#agn galaxy plot 
#plt.scatter(agn_df["z_Spec"], agn_df["O3_Hb_ratio"], color = "red") #plot the redshift vs the OIII/H-beta ratio to see if there is any relationship between them
#plt.yscale("log")
 #set the y-axis to a logarithmic scale bc emission line ratios
# span several orders of magnitude
#plt.xlabel('Redshift')
#plt.ylabel('[OIII] /H-beta Ratio')
#plt.title('AGN Diagnostics in Early Galaxies')
#plt.savefig("figures/o3_hb_ratio_vs_redshift.png", dpi=300, bbox_inches = "tight")
#plt.show()


strong_sf = df.sort_values(
    by= "HA_6563_flux",
    ascending = False
)

#print(sf_df.shape) 
#print(agn_df.shape)

#print(sf_df.head(10))
#ok so this gave the tables full of all the neccesary info on the SF galaxies 
#now for agn 
#print(agn_df.head(10))
#ok fire it works yippee 

# after moving the both datasets that were printed out to google sheets, 
# I picked a few galaxies that looked interesting to analyze 
# 4873, is a great candidate for being both star forming and AGN 
# 2430 is a strong star forming comparison galaxy
# 4297 is a high redshift galaxy that is a strong agn candidate
# 4404 is a strong starforming + good 0II + OIII lines (balanced)


#okay going to attempt to use ASTROQUERY to find the spectral data for these galaxies
#organized galaxie id #'s. their ra and dec, and now a loop to query for the spectra of each galaxy and download it
galaxies= [4873, 2430, 4297, 4404]
ra_list = [53.166109, 53.128186, 53.155788,53.115372 ]
dec_list= [-27.785739, -27.787688, -27.815202, -27.814771]

spectra_files = {}

for i in range(len(galaxies)):
    print(f"Querying for galaxy {galaxies[i]} with RA: {ra_list[i]} and Dec: {dec_list[i]}")

    obs= Observations.query_criteria(obs_collection = "HLSP", 
                                 coordinates = f"{ra_list[i]} {dec_list[i]}",
                                 radius = '2 arcsec')
    obs_pd = obs.to_pandas()

    jades_obs_pd = obs_pd[
    (obs_pd["project"].astype(str).str.contains("JADES", case=False, na=False)) |
    (obs_pd["target_name"].astype(str).str.contains("JADES", case=False, na=False)) |
    (obs_pd["instrument_name"].astype(str).str.contains("NIRSPEC", case=False, na=False))
    ]

    jades_obs = Table.from_pandas(jades_obs_pd)

    products = Observations.get_product_list(jades_obs)

    matches = products[
        [("x1d" in str(name).lower()) or ("spec1d" in str(name).lower()) or ("prism" in str(name).lower())
          for name in products["productFilename"]]
    ]

    print(matches["productFilename"])
    prism_x1d = matches [
        [("clear-prism" in str(name).lower()) and ("x1d" in str(name).lower())
         for name in matches["productFilename"]]]

    goodies = Observations.download_products(prism_x1d)

    local_path = goodies["Local Path"][0]
    print("Downloaded:" , local_path)

    spectra_files[galaxies[i]] = local_path
print(spectra_files)

#going to individually organize their files in order to plot their spectra 
spectra_data = {}
for galaxy_id in spectra_files:
    fname = spectra_files[galaxy_id]

    hdu = fits.open(fname)
    data = hdu[1].data

    spectra_data[galaxy_id] = {
        "wavelength": data["WAVELENGTH"],
        "flux": data["FLUX"],
        "flux_err": data["FLUX_ERR"]
    }

print(spectra_data.keys())
 
# emission lines 
emission_lines = {
    "H\u03B1 6563" : 0.6563,
    "H\u03B2 4861" : 0.4861,
    "[OII] 3727" : 0.3727,
    "[OIII 5007" : 0.5007,
}

# Plot settings
galaxy_colors = {
    4873: "tab:blue",
    2430: "tab:green",
    4297: "tab:red",
    4404: "tab:purple",
}

galaxy_notes = {
    4873: "Strong SF + strong ionization",
    2430: "Star-forming comparison",
    4297: "High-redshift ionized source",
    4404: "SF + OII/OIII detected",
}

fig, axes = plt.subplots(
    nrows=4,
    ncols=1,
    figsize=(11, 10),
    sharex=True
)

#first lets account for redshift to get each spectra in its rest wavelength frame
#lamba_rest = lambda_obs / 1+z 

for index, galaxy_id in enumerate(galaxies):

    z = df[df["NIRSpec_ID"] == galaxy_id]["z_Spec"].values[0]

    wavelength_obs = spectra_data[galaxy_id]["wavelength"]
    wavelength_rest = wavelength_obs / (1 + z)
    flux = spectra_data[galaxy_id]["flux"]

    ax = axes[index]

    # remove noisy far-blue edge
    #didnt work??? -NEED TO FIX 
    mask = wavelength_rest > 0.12

    ax.plot(
        wavelength_rest[mask],
        flux[mask],
        color=galaxy_colors[galaxy_id],
        alpha=0.9
    )

    ax.set_ylabel("Flux")
    ax.set_title(
        f"Galaxy {galaxy_id}: {galaxy_notes[galaxy_id]} (z={z:.2f})"
    )

    # Add emission-line markers to EACH panel
    for label, line_wave in emission_lines.items():
        ax.axvline(
            x=line_wave,
            color="black",
            linestyle="--",
            alpha=0.45
        )

        ax.text(
            line_wave,
            ax.get_ylim()[1] * 0.85,
            label,
            rotation=90,
            fontsize=8,
            verticalalignment="top"
        )

axes[-1].set_xlabel("Rest-frame Wavelength ($\\mu$m)")

fig.suptitle(
    "Rest-frame JADES NIRSpec Prism Spectra of Selected High-Redshift Galaxies",
    y=1.02
)
plt.tight_layout()
plt.savefig("figures/rest_frame_spectra_four_galaxies.png", dpi=300, bbox_inches="tight")
plt.show()

#YAY it worked 
#data anylsis to be continued
