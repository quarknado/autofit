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

xinit, yinit = ad.file_reader(directory + '/' + f)


spe = plt.figure(figsize = (15,7))
ax = spe.add_subplot(111)
ax.plot(xinit, yinit)
ax.set_xticks(np.arange(0,max(xinit), 100))
spe.show()

'''
###########################################Contaminant Deletion##################################################
'''

cont_del = 'y'# ad.y_or_n("Here is your spectrum. Would you like to remove any contaminants?")

if cont_del == 'y':
    xclip, yclip = [xinit, yinit]

    while True:
        
        xclipbuffer, yclipbuffer = [xclip, yclip]
        plt.plot(xclip, yclip)
        plt.show()        
        xclip, yclip = sig.contaminantclipper(xclip, yclip)
        plt.plot(xclip, yclip)
        plt.show()
        confirmclip = ad.y_or_n('Are you happy with this removal?')

        if confirmclip == 'n':

            cancel = ad.y_or_n('Are you sure there are any contaminants?')
            if cancel == 'n':
                xclip, yclip = [xinit, yinit]
                break
            if cancel == 'y':
                xclip, yclip = [xclipbuffer, yclipbuffer]
                continue
            
        else:
            pass

    
        moreclip = ad.y_or_n('Would you like to remove any more contaminants?')
        if moreclip == 'y':
            continue
        else:
            break

plt.plot(xclip,yclip)

















