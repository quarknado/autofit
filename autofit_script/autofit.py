## adminfunctions.py
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

directory = os.getcwd() + '/spectra'

fileno = ad.listfiles(directory)
fileselect = 1#ad.inputfileselector(fileno)
f = ad.openfilefromlist(fileselect, directory)

xinit, yinit = ad.file_reader(directory + '/' + f)

print("Here is your spectrum. Would you like to remove any contaminants?(y/n)")

spe = plt.figure(figsize = (15,7))
ax = spe.add_subplot(111)
ax.plot(xinit, yinit)
ax.set_xticks(np.arange(0,max(xinit), 100))
plt.show()
