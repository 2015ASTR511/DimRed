# -*- coding: utf-8 -*-
"""
Code to read in apStar fits files and plot a subset of combined spectra.
Run EXTREMELY NAIVE PCA on the "learning dataset" of 265 apStar files.
  
Author: Grace Telford
Date: 2/18/2015
For Dimensionality Reduction Project, ASTR 511 
"""

import numpy as np
import pyfits #maybe should use astropy.io instead?? 
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
plt.ion()

# read in metallicities from file from CAS aligned line by line with spectra filenames
# check how many disk vs. halo stars we have; get rid of stars with no metallicity determination

metallicities = np.genfromtxt('metallicity.txt',delimiter=',',usecols=0)
good_z = metallicities != -9999.
metallicities = metallicities[good_z]

condHalo = (metallicities > -2.5) * (metallicities < -1.1)
condDisk = (metallicities > -0.9) * (metallicities < 0.6) 

print 'Number of halo stars: ', np.sum(condHalo)
print 'Number of disk stars: ', np.sum(condDisk)



filenames = np.loadtxt('spec_filenames.txt', dtype=str)
filenames = filenames[good_z]

N_spec = len(filenames)

# just open one apStar file for now

file = pyfits.open('APOGEE/'+filenames[1])

# define wavelength grid - should be same for all apStar files

start = file[0].header['CRVAL1']
delta = file[0].header['CDELT1']
N_wave = file[0].header['NWAVE']
end = start + N_wave * delta

wavelength = np.logspace(start, end, N_wave)

# make simple plot of first spectrum

units = 10**-17 # conversion required to get physical flux units of erg/s/cm^2/Ang
flux = file[1].data[0,:] * units
error = file[2].data[0,:] * units

file.close()

plt.figure()
plt.plot(wavelength, flux)
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel(r'Flux (erg/cm$^2$/s/$\AA$)')

plt.figure()
plt.plot(wavelength, flux, label='Flux')
plt.plot(wavelength, error, label='Error')
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel(r'Flux (erg/cm$^2$/s/$\AA$)')
plt.yscale('log')
plt.legend()

# create array of all 256 spectra

spectra = np.zeros((N_spec, N_wave))

#mask = error < 10**-15.5 # != np.max(error) # mask out places with really big errors/flux=0
#spectra = np.zeros((N_spec, np.sum(mask)))

for index in range(N_spec):
    file = pyfits.open('APOGEE/'+filenames[index])
    if file[0].header['NVISITS'] == 1:
        spectra[index,:] = file[1].data #* units
        #spectra[index,:] = file[1].data[mask]
    else:
        spectra[index,:] = file[1].data[0,:] #* units
        #spectra[index,:] = file[1].data[0,mask]
    file.close()
    
# run naive PCA without masking anything just to see if it works
    
pca = PCA(n_components=5)
pca.fit(spectra)
print 'Variances of recovered PCs: ', pca.explained_variance_ratio_

#wavelength = wavelength[mask]

plt.figure(figsize=(10,8))
plt.title('All Spectra')

plt.subplot(311)
#plt.xlabel('Wavelength ($\mathrm{\AA}$)')
plt.ylabel('Mean')
plt.xlim([15000, 17000])
plt.plot(wavelength, pca.mean_)
    
plt.subplot(312)
#plt.xlabel('Wavelength ($\mathrm{\AA}$)')
plt.ylabel('PC 1')
plt.xlim([15000, 17000])
plt.text(15200, 0.004, 'Variance: %g' % pca.explained_variance_ratio_[0])
plt.plot(wavelength, pca.components_[0])
    
plt.subplot(313)
plt.xlabel('Wavelength ($\mathrm{\AA}$)')
plt.ylabel('PC 2')
plt.xlim([15000, 17000])
plt.text(15200, -0.04, 'Variance: %g' % pca.explained_variance_ratio_[1])
plt.plot(wavelength, pca.components_[1])