"""
Code to read in the allStar .fits table and apStar .fits files containing
combined spectra and save arrays of all spectra and associated errors and
star flags.

Spectra will be saved in chunks spanning 100K in Teff. Also for each chunk
an array of indices in the allStar file (excluding rows where LOCATION_ID = 1)
corresponding to those spectra is saved; these masks will be used to select
relevant parameters and flags later.

Author: Grace Telford
Date: 3/2/15
For Dimensionality Reduction Project, ASTR 511
"""

import numpy as np
import pyfits

allStar = pyfits.open('/astro/store/scratch/tmp/DimRed/allStar-v603.fits')

fields = allStar[1].data['LOCATION_ID']
""" 
~900 objects have field set to 1; wget did not find spectra for these so we
are excluding those stars.
"""
mask = fields != 1
Teff = allStar[1].data['TEFF'][mask]
allStar.close()

filenames = np.loadtxt('/astro/store/scratch/tmp/DimRed/apStar_files.txt', dtype=str)

# open first apStar file to get number of wavelength samples
file = pyfits.open('/astro/store/scratch/tmp/DimRed/apStar/'+filenames[0])
N_wave = file[0].header['NWAVE']
file.close() 

# loop over Teff bins and save arrays of spectra, flags, and masks for each

T_bins = np.arange(3600, 5600, 100)

for T in T_bins:

    temp_inds = (Teff > T) * (Teff < (T + 100.))
    N_spec = np.sum(temp_inds)

    file_chunk = filenames[temp_inds]

    # generate array of spectra and errors
    spectra = np.zeros((2, N_spec, N_wave))
    flags = np.zeros((N_spec))

    for index in range(N_spec):
        file = pyfits.open('/astro/store/scratch/tmp/DimRed/apStar/'+file_chunk[index])

        # save bitwise or star flag for all visits to each star
        flags[index] = file[0].header['STARFLAG']

        if file[0].header['NVISITS'] == 1:
            spectra[0, index, :] = file[1].data #* units
            spectra[1, index, :] = file[2].data

        else:
            spectra[0, index, :] = file[1].data[0, :] #* units
            spectra[1, index, :] = file[2].data[0, :]

        file.close()

    # save all arrays as .npz files
    np.save('mask'+str(T), temp_inds) 
    np.save('spectra'+str(T), spectra)
    np.save('flags'+str(T), flags)
