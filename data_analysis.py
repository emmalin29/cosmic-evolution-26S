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
    'O3_4959_flux', #flux of the OIII 4959 line -AGN indicator
    'O3_4959_err', #error on the OIII 4959 flux measurement

]
#future use for JADES viewer pics (potentially)
print(df['NIRCam_ID'])

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

#——————————————————————————————————————————#
#phase 2: 
#going to use ASTROQUERY to find the spectral data for these galaxies
#organized galaxie id #'s. their ra and dec, and now a loop to query for the spectra of each galaxy and download it
galaxies= [4873, 2430, 4297, 4404]
ra_list = [53.166109, 53.128186, 53.155788,53.115372 ]
dec_list= [-27.785739, -27.787688, -27.815202, -27.814771]


galaxies = [4873, 2430, 4297, 4404]

selected_info = df[df["NIRSpec_ID"].isin(galaxies)][
    [
        "NIRSpec_ID",
        "NIRCam_ID",
        "RA_TARG",
        "Dec_TARG",
        "z_Spec",
        "z_Spec_flag"
    ]
]

print(selected_info)

spectra_files = {}

"""for i in range(len(galaxies)):
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
print(spectra_files)"""
spectra_files = {
    4873: "./mastDownload/HLSP/hlsp_jades_jwst_nirspec_goods-s-mediumhst-00004873_clear-prism_v1.0/hlsp_jades_jwst_nirspec_goods-s-mediumhst-00004873_clear-prism_v1.0_x1d.fits" ,
    2430: "./mastDownload/HLSP/hlsp_jades_jwst_nirspec_goods-s-mediumhst-00002430_clear-prism_v1.0/hlsp_jades_jwst_nirspec_goods-s-mediumhst-00002430_clear-prism_v1.0_x1d.fits",
    4297: "./mastDownload/HLSP/hlsp_jades_jwst_nirspec_goods-s-deephst-00004297_clear-prism_v1.0/hlsp_jades_jwst_nirspec_goods-s-deephst-00004297_clear-prism_v1.0_x1d.fits",
    4404: "./mastDownload/HLSP/hlsp_jades_jwst_nirspec_goods-s-deephst-00004404_clear-prism_v1.0/hlsp_jades_jwst_nirspec_goods-s-deephst-00004404_clear-prism_v1.0_x1d.fits"
}


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

#print(spectra_data.keys())
 
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
    figsize=(15, 8),
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
    mask = wavelength_rest > 0.2

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
            ax.get_ylim()[1] * 0.85+ 0.003,
            label,
            rotation=90,
            fontsize=9,
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

#——————————————————————————————————————————#
#Formal spectra line analysis to see if we got an agn + SF candidates 

#—————————#
#only going to look at Hb, OIII, and HA, since OII's peaks aren't very clear 
#—————————#

plt.rcParams['figure.figsize'] = (12.0, 8.0)

#4873 first 
#Load up  spectrum
glx4873_spectrum = spectra_data[4873]
flux = glx4873_spectrum['flux']
wavelength = glx4873_spectrum['wavelength']

pixels_full = wavelength / (1 + df[df["NIRSpec_ID"] == 4873]["z_Spec"].values[0]) #convert to rest frame wavelength by dividing by 1+z
signal_full = flux

#plt.plot(pixels_full,signal_full)  #fill in what to plot, then run this cell again to see if it matches
x1= 0.45
x2= 0.7
mask = (pixels_full >= x1) & (pixels_full <= x2)
pixels = pixels_full[mask]
signal = signal_full[mask]

"""
plt.xlim(x1,x2) 
plt.ylim(-0.2e-19, 1.2e-19)
plt.xlabel(' wavelength (microns)')
plt.ylabel('flux')
plt.show()
"""
#Find the peaks 

threshold =1.5e-20 #You can just pick slightly lower than the lowest peak you want to centroid
peaks = []    #x positions of the peaks, or rather, their index
for i in range(2,len(signal)-2): #len(signal)-1 because you will be checking the value after than your last i 
    if (signal[i] > signal[i-1]) and (signal[i] > signal[i+1]) and (signal[i] > threshold) :  #three conditions to be a peak
        if((signal[i] > signal[i-2]) and (signal[i] > signal[i+2])):
            peaks.append(i)

c= 299792.458 #speed of light in km/s 

centroids = [] #Values for all the centroids
xminlist = [] 
xmaxlist = [] 

print("Galaxy 4873 FWHM DATA: ")
for i in range(len(peaks)):
    #Calculate how far backward and forward to go:
    half_max = signal[peaks[i]]/2
    indices= np.arange(len(signal))
    xmin =  int(indices[(indices<=peaks[i])&(signal<=half_max)][-1])
    xmax =  int(indices[(indices>=peaks[i])&(signal<=half_max)][0])
    x_range = pixels[xmin:xmax +1 ]
    I_range = signal[xmin:xmax +1 ]
    x_range = np.array(x_range)
    I_range = np.array(I_range)
    

    #finding FWHM of each peak
    lambda_0 = np.sum(x_range * I_range) / np.sum(I_range)
    centroids.append(lambda_0)
    delta_lambda = pixels[xmax] - pixels[xmin] 
    velocity= c * (delta_lambda / lambda_0 )
    xminlist.append(pixels[xmin])
    xmaxlist.append(pixels[xmax]) 

    print (f"Peak {i+1} at {lambda_0:.4f} microns:")
    print(f"FWHM : {delta_lambda:.5f}microns")
    print(f"Velocity: {velocity:.1f} km/s\n")

    
def plot_vert(x): 
    '''
    Just plots vertical lines, in blue dashes
    '''
    plt.axvline(x, color='DarkBlue', ls='-.')


#for i in centroids[1:]: #Call my plotting function on every centroid except the first
#    plot_vert(i)
#plt.axvline(centroids[0], color='DarkBlue', ls='-.', label='Centroid', lw=2) #Reserve the first so I don't have a million "centroid" labels


plt.plot(pixels, signal, 'blue', label='4873 spectrum') #Plot the actual spectrum
ymin, ymax = plt.gca().get_ylim()
plt.axvline(0.4861, ls= "--", color= "lightskyblue", lw = 1)
plt.text( 0.4861 + 0.004, ymax*0.85, "H\u03B2 4861", rotation=90, color="lightskyblue", fontsize=11, va="top")
plt.axvline(0.5007,  ls= "--", color= "lightseagreen", lw= 1)
plt.text( 0.5007 + 0.004, ymax*0.85, "[OIII] 5007", rotation=90, color="lightseagreen", fontsize=11, va="top")
plt.axvline(0.6563, ls= "--", color= "red", lw = 1)
plt.text( 0.6563 + 0.004, ymax*0.85, "H\u03B1 6563", rotation=90, color="red", fontsize=11, va="top")

for i, xmin in enumerate(xminlist):
    if i == 0: 
        plt.axvline(xmin, color= "grey", ls=":", alpha = 0.8, lw = 1.8 , label = "FWHM Boundaries") 
    else:
        plt.axvline(xmin, color= "grey", ls=":", alpha = 0.8, lw = 1.8) 
for xmax in xmaxlist:
    plt.axvline(xmax, color= "grey", ls=":", alpha = 0.8, lw =1.8 ) 

plt.title("FWHM of Galaxy 4873: AGN + SF Candidate")
plt.xlabel('Rest-frame Wavelength ($\\mu$m)')
plt.ylabel('Flux')
plt.legend(loc=2)
plt.show()
print(centroids)
#——————————————————————————————————————————#
#galaxy 2430- SF
glx2430_spectrum = spectra_data[2430]
flux = glx2430_spectrum['flux']
wavelength = glx2430_spectrum['wavelength']
pixels_full = wavelength / (1 + df[df["NIRSpec_ID"] == 2430]["z_Spec"].values[0]) 
signal_full = flux

#plt.plot(pixels,signal) 
x1= 0.45 
x2= 0.7 
mask = (pixels_full >= x1) & (pixels_full <= x2)
pixels = pixels_full[mask]
signal = signal_full[mask]

#Find the peaks 
threshold =1.0e-20 
peaksSF = []  
for i in range(2,len(signal)-2): 
    if (signal[i] > signal[i-1]) and (signal[i] > signal[i+1]) and (signal[i] > threshold) :  
        if((signal[i] > signal[i-2]) and (signal[i] > signal[i+2])):
            peaksSF.append(i)

c= 299792.458 

centroidsSF = [] 
xminlistSF = [] 
xmaxlistSF = [] 

print("Galaxy 2430 (SF) FWHM DATA: ")
for i in range(len(peaksSF)):
    half_max = signal[peaksSF[i]]/2
    indices= np.arange(len(signal))
    xmin =  int(indices[(indices<=peaksSF[i])&(signal<=half_max)][-1])
    xmax =  int(indices[(indices>=peaksSF[i])&(signal<=half_max)][0])
    x_range = pixels[xmin:xmax +1 ]
    I_range = signal[xmin:xmax +1 ]
    x_range = np.array(x_range)
    I_range = np.array(I_range)


    #finding FWHM of each peak
    lambda_0 = np.sum(x_range * I_range) / np.sum(I_range)
    centroidsSF.append(lambda_0)
    delta_lambda = pixels[xmax] - pixels[xmin] 
    velocity= c * (delta_lambda / lambda_0 )
    xminlistSF.append(pixels[xmin])
    xmaxlistSF.append(pixels[xmax]) 

    print (f"Peak {i+1} at {lambda_0:.4f} microns:")
    print(f"FWHM : {delta_lambda:.5f}microns")
    print(f"Velocity: {velocity:.1f} km/s\n")

    
#for i in centroidsSF[1:]: 
#    plot_vert(i)
#plt.axvline(centroidsSF[0], color='DarkBlue', ls='-.', label='Centroid')

plt.plot(pixels, signal, 'green', label='Spectrum') 
ymin, ymax = plt.gca().get_ylim()
plt.axvline(0.4861, ls= "--", color= "lightskyblue", lw = 1)
plt.text( 0.4861 + 0.004, ymax*0.85, "H\u03B2 4861", rotation=90, color="lightskyblue", fontsize=11, va="top")
plt.axvline(0.5007,  ls= "--", color= "lightseagreen", lw= 1)
plt.text( 0.5007 + 0.004, ymax*0.85, "[OIII] 5007", rotation=90, color="lightseagreen", fontsize=11, va="top")
plt.axvline(0.6563, ls= "--", color= "red", lw = 1)
plt.text( 0.6563 + 0.004, ymax*0.85, "H\u03B1 6563", rotation=90, color="red", fontsize=11, va="top")

for i, xmin in enumerate(xminlistSF):
    if i == 0: 
        plt.axvline(xmin, color= "grey", ls=":", alpha = 0.8, lw = 1.8 , label = "FWHM Boundaries") 
    else:
        plt.axvline(xmin, color= "grey", ls=":", alpha = 0.8, lw = 1.8) 
for xmax in xmaxlistSF:
    plt.axvline(xmax, color= "grey", ls=":", alpha = 0.8, lw =1.8 ) 
plt.title("FWHM of Galaxy 2430: SF")
plt.xlabel('Rest-frame Wavelength ($\\mu$m)')
plt.ylabel('Flux')
plt.legend(loc=2)
plt.show()
print(centroidsSF)

#——————————————————————————————————————————#
#galaxy 4297- AGN 
glx4297_spectrum = spectra_data[4297]
flux = glx4297_spectrum['flux']
wavelength = glx4297_spectrum['wavelength']
pixels_full = wavelength / (1 + df[df["NIRSpec_ID"] == 4297]["z_Spec"].values[0]) #convert to rest frame wavelength by dividing by 1+z
signal_full = flux

#plt.plot(pixels,signal) #fill in what to plot, then run this cell again to see if it matches
x1= 0.45 
x2= 0.7 
mask = (pixels_full >= x1) & (pixels_full <= x2)
pixels = pixels_full[mask]
signal = signal_full[mask]

#Find the peaks 
threshold =0.65e-20 #You can just pick slightly lower than the lowest peak you want to centroid
peaksAGN = []    #x positions of the peaks, or rather, their index
for i in range(2,len(signal)-2): #len(signal)-1 because you will be checking the value after than your last i 
    if (signal[i] > signal[i-1]) and (signal[i] > signal[i+1]) and (signal[i] > threshold) :  #three conditions to be a peak
        if((signal[i] > signal[i-2]) and (signal[i] > signal[i+2])):
            peaksAGN.append(i)


c= 299792.458 #speed of light in km/s 

centroidsAGN = [] #Values for all the centroids
xminlistAGN = [] 
xmaxlistAGN = [] 

print("Galaxy 4297(AGN) FWHM DATA: ")
for i in range(len(peaksAGN)):
    #Calculate how far backward and forward to go:
    half_max = signal[peaksAGN[i]]/2
    indices= np.arange(len(signal))
    xmin =  int(indices[(indices<=peaksAGN[i])&(signal<=half_max)][-1])
    xmax =  int(indices[(indices>=peaksAGN[i])&(signal<=half_max)][0])
    x_range = pixels[xmin:xmax +1 ]
    I_range = signal[xmin:xmax +1 ]
    x_range = np.array(x_range)
    I_range = np.array(I_range)


    #finding FWHM of each peak
    lambda_0 = np.sum(x_range * I_range) / np.sum(I_range)
    centroidsAGN.append(lambda_0)
    delta_lambda = pixels[xmax] - pixels[xmin] 
    velocity= c * (delta_lambda / lambda_0 )
    xminlistAGN.append(pixels[xmin])
    xmaxlistAGN.append(pixels[xmax]) 

    print (f"Peak {i+1} at {lambda_0:.4f} microns:")
    print(f"FWHM : {delta_lambda:.5f}microns")
    print(f"Velocity: {velocity:.1f} km/s\n")

     
#for i in centroidsAGN[1:]: 
#    plot_vert(i)
#plt.axvline(centroidsAGN[0], color='DarkBlue', ls='-.', label='Centroid')

plt.plot(pixels, signal, 'red', label='Spectrum') 
ymin, ymax = plt.gca().get_ylim()
plt.axvline(0.4861, ls= "--", color= "lightskyblue", lw = 1)
plt.text( 0.4861 + 0.004, ymax*0.85, "H\u03B2 4861", rotation=90, color="lightskyblue", fontsize=11, va="top")
plt.axvline(0.5007,  ls= "--", color= "lightseagreen", lw= 1)
plt.text( 0.5007 + 0.004, ymax*0.85, "[OIII] 5007", rotation=90, color="lightseagreen", fontsize=11, va="top")
plt.axvline(0.6563, ls= "--", color= "darkred", lw = 1)
plt.text( 0.6563 + 0.004, ymax*0.85, "H\u03B1 6563", rotation=90, color="darkred", fontsize=11, va="top")

for i, xmin in enumerate(xminlistAGN):
    if i == 0: 
        plt.axvline(xmin, color= "grey", ls=":", alpha = 0.8, lw = 1.8 , label = "FWHM Boundaries") 
    else:
        plt.axvline(xmin, color= "grey", ls=":", alpha = 0.8, lw = 1.8) 
for xmax in xmaxlistAGN:
    plt.axvline(xmax, color= "grey", ls=":", alpha = 0.8, lw =1.8 ) 
plt.title("FWHM of Galaxy 4297: AGN")
plt.xlabel('Rest-frame Wavelength ($\\mu$m)')
plt.ylabel('Flux')
plt.legend(loc=2)
plt.show()
print(centroidsAGN)

#——————————————————————————————————————————#
#galaxy 4404- Balanced  
glx4404_spectrum = spectra_data[4404]
flux = glx4404_spectrum['flux']
wavelength = glx4404_spectrum['wavelength']
pixels_full = wavelength / (1 + df[df["NIRSpec_ID"] == 4404]["z_Spec"].values[0]) #convert to rest frame wavelength by dividing by 1+z
signal_full = flux

#plt.plot(pixels,signal) #fill in what to plot, then run this cell again to see if it matches
x1= 0.45 
x2= 0.7 
mask = (pixels_full >= x1) & (pixels_full <= x2)
pixels = pixels_full[mask]
signal = signal_full[mask]

#Find the peaks 
threshold =1.0e-20 #You can just pick slightly lower than the lowest peak you want to centroid
peaksBL = []    #x positions of the peaks, or rather, their index
for i in range(2,len(signal)-2): #len(signal)-1 because you will be checking the value after than your last i 
    if (signal[i] > signal[i-1]) and (signal[i] > signal[i+1]) and (signal[i] > threshold) :  #three conditions to be a peak
        if((signal[i] > signal[i-2]) and (signal[i] > signal[i+2])):
            peaksBL.append(i)


c= 299792.458 #speed of light in km/s 

centroidsBL = [] #Values for all the centroids
xminlistBL = [] 
xmaxlistBL = [] 

print("Galaxy 4404 (Balanced) FWHM DATA: ")
for i in range(len(peaksBL)):
    #Calculate how far backward and forward to go:
    half_max = signal[peaksBL[i]]/2
    indices= np.arange(len(signal))
    xmin =  int(indices[(indices<=peaksBL[i])&(signal<=half_max)][-1])
    xmax =  int(indices[(indices>=peaksBL[i])&(signal<=half_max)][0])
    x_range = pixels[xmin:xmax +1 ]
    I_range = signal[xmin:xmax +1 ]
    x_range = np.array(x_range)
    I_range = np.array(I_range)


   #finding FWHM of each peak
    lambda_0 = np.sum(x_range * I_range) / np.sum(I_range)
    centroidsBL.append(lambda_0)
    delta_lambda = pixels[xmax] - pixels[xmin] 
    velocity= c * (delta_lambda / lambda_0 )
    xminlistBL.append(pixels[xmin])
    xmaxlistBL.append(pixels[xmax]) 

    print (f"Peak {i+1} at {lambda_0:.4f} microns:")
    print(f"FWHM : {delta_lambda:.5f}microns")
    print(f"Velocity: {velocity:.1f} km/s\n")

    
#for i in centroidsBL[1:]: 
#    plot_vert(i)
#plt.axvline(centroidsBL[0], color='DarkBlue', ls='-.', label='Centroid')

plt.plot(pixels, signal, 'purple', label='Spectrum') 
ymin, ymax = plt.gca().get_ylim()
plt.axvline(0.4861, ls= "--", color= "lightskyblue", lw = 1)
plt.text( 0.4861 + 0.004, ymax*0.85, "H\u03B2 4861", rotation=90, color="lightskyblue", fontsize=11, va="top")
plt.axvline(0.5007,  ls= "--", color= "lightseagreen", lw= 1)
plt.text( 0.5007 + 0.004, ymax*0.85, "[OIII] 5007", rotation=90, color="lightseagreen", fontsize=11, va="top")
plt.axvline(0.6563, ls= "--", color= "red", lw = 1)
plt.text( 0.6563 + 0.004, ymax*0.85, "H\u03B1 6563", rotation=90, color="red", fontsize=11, va="top")

for i, xmin in enumerate(xminlistBL):
    if i == 0: 
        plt.axvline(xmin, color= "grey", ls=":", alpha = 0.8, lw = 1.8 , label = "FWHM Boundaries") 
    else:
        plt.axvline(xmin, color= "grey", ls=":", alpha = 0.8, lw = 1.8) 
for xmax in xmaxlistBL:
    plt.axvline(xmax, color= "grey", ls=":", alpha = 0.8, lw =1.8 ) 
plt.title("FWHM of Galaxy 4404: Balanced ")
plt.xlabel('Rest-frame Wavelength ($\\mu$m)')
plt.ylabel('Flux')
plt.legend(loc=2)
plt.show()
print(centroidsBL)

#——————————————————————————————————————————#
