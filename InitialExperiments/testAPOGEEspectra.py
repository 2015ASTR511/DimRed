# -*- coding: utf-8 -*-
"""
Code to read in apStar fits files and plot a subset of combined spectra.
Run PCA on the "learning dataset" of 265 apStar files.
  
Author: Grace Telford, Diana Windemuth
Date: 2/18/2015
For Dimensionality Reduction Project, ASTR 511 

Updated 2/22/2015 by Diana Windemuth
Updated 2/25/2015 by Grace Telford
"""

import numpy as np
import pyfits #maybe should use astropy.io instead?? 
from sklearn.decomposition import PCA
import matplotlib as mpl
mpl.rc('xtick', labelsize=14)     
mpl.rc('ytick', labelsize=14)  
mpl.rc('axes', labelsize=16)
mpl.rc('legend', fontsize=16)
#mpl.rc('lines', linewidth=2)
mpl.rc('font', size=16, family='serif', style='normal', variant='normal', stretch='normal', weight='normal')
import matplotlib.pyplot as plt
plt.ion()

# read in metallicities from file from CAS aligned line by line with spectra filenames
# check how many disk vs. halo stars we have; get rid of stars with no metallicity determination

def readin(metname = 'metallicity.txt', specname = 'spec_filenames.txt', plot=True):
    metallicities = np.genfromtxt(metname,delimiter=',',usecols=0)
    good_z = metallicities != -9999.
    metallicities = metallicities[good_z]
    
    condHalo = (metallicities > -2.5) * (metallicities < -1.1)
    condDisk = (metallicities > -0.9) * (metallicities < 0.6) 
    
    print 'Number of halo stars: ', np.sum(condHalo)
    print 'Number of disk stars: ', np.sum(condDisk)
    
    filenames = np.loadtxt(specname, dtype=str)
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
    
    units = 1e-17 # conversion required to get physical flux units of erg/s/cm^2/Ang
    flux = file[1].data[0,:] * units
    error = file[2].data[0,:] * units
    
    file.close()
    
    if plot:
        plt.figure()
        plt.plot(wavelength, flux, label='Flux')
        plt.plot(wavelength, error, label='Error')
        plt.xlabel(r'Wavelength ($\AA$)')
        plt.ylabel(r'Flux (erg/cm$^2$/s/$\AA$)')
        plt.yscale('log')
        plt.legend()

    
    # create array of all 256 spectra
    spectra = np.zeros((2, N_spec, N_wave))    
    flags = np.zeros((N_spec))

    
    for index in range(N_spec):
        file = pyfits.open('APOGEE/'+filenames[index])
        
        # save bitwise or star flag for all visits to each star
        flags[index] = file[0].header['STARFLAG']
        
        if file[0].header['NVISITS'] == 1:
            spectra[0, index, :] = file[1].data #* units
            spectra[1, index, :] = file[2].data

        else:
            spectra[0, index, :] = file[1].data[0, :] #* units
            spectra[1, index, :] = file[2].data[0, :]

        file.close()
        
    # remove spectra where starflag was nonzero for any visit; can refine later. Binary digits 0, 3, 4 == BAD
    spectra = spectra[:, flags == 0, :]
    print 'Number of spectra after removing those with flags set: ', np.shape(spectra)[1]
    
    return wavelength, spectra

def naivePCA(wavelength, spectra, mask=None, n_components=5, plot=True, figtitle=None):

    pca = PCA(n_components=n_components)

    if mask is not None:
        spectra = spectra[:, :, mask]
        wavelength = wavelength[mask]

    pca.fit(spectra[0, :, :])
    print 'Variances of recovered PCs: ', pca.explained_variance_ratio_

    if plot:
        plt.figure(figsize=(12,8))
        plt.suptitle(figtitle)
        
        plt.subplot(3,2,1)
        plt.ylabel('Mean')
        plt.xlim([15000, 17000])
        plt.plot(wavelength, pca.mean_)

        for ii in range(n_components):
            plt.subplot(3,2,ii+2)
            plt.ylabel('PC '+str(ii+1))
            plt.xlim([15000, 17000])
            plt.text(15200, np.min(pca.components_[ii]), 'Variance: %g' % pca.explained_variance_ratio_[ii],fontsize=12)
            plt.plot(wavelength, pca.components_[ii])
        #plt.tight_layout()

        # make scree plot
        plt.figure(figsize=(8,6))
        plt.plot(np.arange(1,6), pca.explained_variance_ratio_)
        plt.title(figtitle)
        plt.xlabel('PC')
        plt.ylabel('Variance')
            

def testMasks(wavelength, spectra):
    # identify common chip gaps and fill them in w/a value 
    # drawn from a gaussian distribution
    meanerr = np.mean(spectra[1, :, :], axis=0)
    gaps = (meanerr == 1e10)
    spectra2 = spectra.copy()
    spectra2[0, :, gaps] = np.median(spectra[0, :, ~gaps], axis=0) * \
                                np.random.normal(1, 0.1, int(np.sum(gaps)))[:, np.newaxis]
    
    naivePCA(wavelength, spectra2, figtitle='Gaussian Gap Filling')
    
    # identify all bad data and fill them with value 
    # drawn from normal distribution
    spectra2 = spectra.copy()
    mask = np.zeros(spectra2.shape[1:], dtype=bool)
    for ii in range(spectra2.shape[1]):
        mask[ii] = (spectra[1, ii, :] == 1e10)
        spectra2[0, ii, mask[ii]] = np.median(spectra[0, ii, ~mask[ii]]) * \
                                np.random.normal(1, 0.1, int(np.sum(mask[ii])))
    
    naivePCA(wavelength, spectra2, figtitle='Gaussian Bad Data Replacement')
    
    # identify common mask where bad data exists for more than 5% of population
    mask = (spectra[1, :, :] == 1e10)
    finalmask = (np.sum(mask, axis=0) > spectra.shape[1]*.05)
    
    naivePCA(wavelength, spectra, ~finalmask, figtitle='Bad Data Masked')
    
    # create mask of chip gaps using wavelengths given on the "using apogee spectra" webpage
    #chipgaps = (wavelength < 15140.) + ((wavelength > 15810.)*(wavelength < 15850.)) \
    #            + ((wavelength > 16440.)*(wavelength < 16470.)) + (wavelength > 16960.)
    
    # now pad the chip gaps by 10 angstroms on each side to make sure big changes in spectrum are removed
    chipgaps = (wavelength < 15170.) + ((wavelength > 15800.)*(wavelength < 15860.)) \
                + ((wavelength > 16430.)*(wavelength < 16480.)) + (wavelength > 16950.)
    
    naivePCA(wavelength, spectra, ~chipgaps, figtitle='Chip Gaps Masked')  


# apply PCA to raw data
wavelength, spectra = readin(plot=False)
naivePCA(wavelength, spectra, figtitle='No Mask')

# apply PCA to masked data
mask = (spectra[1, :, :] == 1e10)
finalmask = (np.sum(mask, axis=0) > spectra.shape[1]*.05)

naivePCA(wavelength, spectra, ~finalmask, figtitle='Bad Data Masked')