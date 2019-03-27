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
    yclip = yinit

    while True:
        
        yclipbuffer = yclip
        
        clip_spectrum = ad.spectrum_plotter(x, yclip)
        clip_spectrum.show()        
        
        yclip = sig.contaminantclipper(x, yclip)
        
        clip_spectrum = ad.spectrum_plotter(x, yclip)
        clip_spectrum.show() 
        
        confirmclip = ad.y_or_n('Are you happy with this removal?')

        if confirmclip == 'n':

            cancel = ad.y_or_n('Are you sure there are any contaminants?')
            if cancel == 'n':
                yclip = yinit
                break
            if cancel == 'y':
                yclip = yclipbuffer
                continue
            
        else:
            pass

    
        moreclip = ad.y_or_n('Would you like to remove any more contaminants?')
        if moreclip == 'y':
            continue
        else:
            break
'''
##########################################Threshold set and fitting###########################################################
'''

width = float(input('What is the approximate FWHM of your peaks in channels'))/(2 * np.sqrt(2*np.log(2)))

ysmooth = sig.smoothe(yclip, width)

smooth_spectrum = ad.spectrum_plotter(x, ysmooth)
smooth_spectrum.show()

thresh = float(input('what threshold would you like to constitute a peak?'))

[xthresh, ythresh] = sig.clipspectrum(x, yclip, ysmooth, thresh)
peak_regions = sig.split_spectrum(xthresh, ythresh, width)
sig.plotall(peak_regions, x, yclip)













