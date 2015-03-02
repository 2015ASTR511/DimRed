"""
Code to read in the allStar .fits table and apStar .fits files containing combined spectra and save arrays of all spectra and relevant parameters (Teff, log(g), [Fe/H], [alpha/M]. Can make other parameter arrays later, just be sure to exlude rows where LOCATION_ID = 1 (no spectra downloaded for those records). No masking is done here, but the bitwise or starflag is saved for each spectrum.

Run in the directory containing the allStar-v603.fits file, or edit the path in the pyfits.open command

Author: Grace Telford
Date: 3/1/15
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
logg = allStar[1].data['LOGG'][mask]
FeH = allStar[1].data['PARAM_M_H'][mask]
alphaM = allStar[1].data['PARAM_ALPHA_M'][mask]

allStar.close()

filenames = np.loadtxt('apStar_files.txt', dtype=str)
N_spec = np.sum(mask)

# open first apStar file to get number of wavelength samples
file = pyfits.open('apStar/'+filenames[0])
N_wave = file[0].header['NWAVE']
file.close()

# generate array of spectra and errors
spectra = np.zeros((2, N_spec, N_wave))
flags = np.zeros((N_spec))

for index in range(N_spec):
    file = pyfits.open('/astro/store/scratch/tmp/DimRed/apStar/'+filenames[index])

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
np.save('Teff', Teff)
np.save('logg', logg)
np.save('FeH', FeH)
np.save('alphaM', alphaM)
np.save('spectra', spectra)
np.save('flags', flags)
