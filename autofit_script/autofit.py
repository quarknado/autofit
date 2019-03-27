## autofit.py
## Automatic fitting routine for transfer reaction spectra
## Ben Cropper 2019

import scipy as sp
import scipy.optimize as op
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
from scipy.signal import find_peaks_cwt
import glob
import os

import adminfunctions as ad
import fittingfunctions as fit
import autosignals as sig


'''
####################################### Read File and Show Spectrum #############################################
'''

directory = os.getcwd() + '/spectra'

fileno = ad.listfiles(directory)
fileselect = 1#ad.inputfileselector(fileno)
f = ad.openfilefromlist(fileselect, directory)

x, yinit = ad.file_reader(directory + '/' + f)


spe = ad.spectrum_plotter(x, yinit)
spe.show()

'''
###########################################Contaminant Deletion##################################################
'''

cont_del = 'y'# ad.y_or_n("Here is your spectrum. Would you like to remove any contaminants?")

if cont_del == 'y':
    yclip = sig.contaminant_clipper(x,yinit)

'''
##########################################Threshold set and Peak Region Detection###########################################################
'''
while True:
    width = float(input('What is the approximate FWHM of your peaks in channels'))/(2 * np.sqrt(2*np.log(2)))

    ysmooth = sig.smoothe(yclip, width)

    smooth_spectrum = ad.spectrum_plotter(x, ysmooth)
    smooth_spectrum.show()

    thresh = float(input('what threshold would you like to constitute a peak?'))

    [xthresh, ythresh] = sig.clipspectrum(x, yclip, ysmooth, thresh)
    peak_regions = sig.split_spectrum(xthresh, ythresh, width)
    sig.plotall(peak_regions, x, yclip)


    recontaminant = ad.y_or_n('Would you like to re-clip for any more contaminants?')
    
    if recontaminant == 'y':
        yclip = sig.contaminant_clipper(x, yclip)
        continue
    else:
        pass

    final_regions = ad.y_or_n('Would you like to re-tune the widths and thresholds for improved peak region detection?')

    if final_regions == 'y':
        continue
    else:
        break


'''
###########################################Fitting##################################
'''









