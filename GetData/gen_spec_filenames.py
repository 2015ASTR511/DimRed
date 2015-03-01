"""
Code to read in the allStar .fits table and generate filenames to feed to the wget command to download all apStar .fits files containing combined spectra.

Run in the directory containing the allStar-v603.fits file, or edit the path in the pyfits.open command

Author: Grace Telford
Date: 2/28/15
For Dimensionality Reduction Project, ASTR 511
"""
import numpy as np
import pyfits

allStar = pyfits.open('allStar-v603.fits')

fields = allStar[1].data['LOCATION_ID']
files = allStar[1].data['FILE']

filenames = []

for index in range(len(fields)):
    filenames.append('http://data.sdss3.org/sas/dr12/apogee/spectro/redux/r5/stars/apo25m/'+str(fields[index])+'/'+files[index])

allStar.close()

filenames = np.array(filenames)
np.savetxt('apStar_filenames.txt', filenames, fmt='%s')
