## autofit.py
## Automatic fitting routine for transfer reaction spectra
## Ben Cropper 2019

import scipy as sp
import scipy.optimize as op
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
from scipy.signal import find_peaks_cwt, find_peaks
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
fileselect = ad.inputfileselector(fileno)
f = ad.openfilefromlist(fileselect, directory)

x, yinit = ad.file_reader(directory + '/' + f)


spe = ad.spectrum_plotter(x, yinit)
spe.show()

'''
###########################################Baseline Subtraction##################################################
'''

base_flag = ad.y_or_n("Does this spectrum need baseline subtraction to fid the peaks.\nThe subtraction will not be applied to the fit. Instead, there will be a linear background applied to each fitting region.")

if base_flag == 'y':
    ybase = sig.baseline_subtraction(x,yinit)
else:
    ybase = yinit


'''
###########################################Contaminant Deletion##################################################
'''

cont_del =  ad.y_or_n("Would you like to remove any contaminants?")

if cont_del == 'y':
    yclipfit, yclipfind  = sig.contaminant_clipper(x,yinit,ybase)
else:
    yclipfit, yclipfind = yinit, ybase

'''
##########################################Threshold set and Peak Region Detection###########################################################
'''
while True:
    spe2 = ad.spectrum_plotter(x, yclipfind, 5)
    spe2.show()
    FWHM = float(input('What is the approximate FWHM of your peaks in channels'))
    width = FWHM/(2 * np.sqrt(2*np.log(2)))


    ysmooth = sig.smoothe(yclipfind, width)

    smooth_spectrum = ad.spectrum_plotter(x, ysmooth)
    smooth_spectrum.show()

    thresh = float(input('what threshold would you like to constitute a peak?'))

    [xthresh, ythresh] = sig.clipspectrum(x, yclipfit, ysmooth, thresh)
    peak_regions = sig.split_spectrum(xthresh, ythresh, width)
    sig.plotall(peak_regions, x, yclipfit)


    recontaminant = ad.y_or_n('Would you like to re-clip for any more contaminants?')
    
    if recontaminant == 'y':
        yclipfit, yclipfind  = sig.contaminant_clipper(x,yclipfit, yclipfind)
        continue
    else:
        pass

    final_regions = ad.y_or_n('Would you like to re-tune the widths and thresholds for improved peak region detection?')

    if final_regions == 'y':
        continue
    else:
        break

fig = plt.figure(figsize = (15,7))
axis = fig.add_subplot(111)
axis.set_xticks(np.arange(0,max(x), 100))

axis.plot(x,yinit, color = 'xkcd:light grey')
for x2,y in peak_regions:
    axis.plot(x2,y)


'''
###########################################Fitting##################################
'''

peak_positions, region_positions = sig.peak_finder(peak_regions, width, FWHM)



pos = np.zeros(len(x))
pos[np.intersect1d(x, peak_positions, return_indices = True)[1]] = 1000
    
axis.plot(x,pos, linewidth = 1)

plt.show()

template, template_pos = ad.template_finder(peak_regions, peak_positions, region_positions,fig)


template_fit = fit.fit(template[0], template[1], template_pos, width, FWHM, rbfix = False, background = True)
ad.printfit(template_fit, template[0], template[1])

rfix, betafix = template_fit[1][-2],template_fit[1][-1]
fitlist = []
for i, region in enumerate(peak_regions):
    xreg = region[0]
    yreg = region[1]    
    ft = fit.fit(xreg, yreg, region_positions[i], width, FWHM, r = rfix, beta = betafix, background = True)
    ad.printfit(ft, xreg, yreg)

    fitlist.append(ft)

'''
##########################################Plot Fits and Save######################
'''

ysub = yclipfit

fitplot = plt.figure(figsize = (15,7))
a = fitplot.add_subplot(111)
a.plot(x, yinit, color = 'xkcd:grey')
a.plot(x, pos, linewidth = 1
)
ofilename = str(f[:-4]) + '_fit.txt'
file = open(ofilename, 'w+')
topline = 'POSITION sPOSITION AREA sAREA WIDTH sWIDTH R BETA\n'
file.write(topline)

for fit in fitlist:
    #plot whole fit    
    a.plot(fit[-2], fit[0], linewidth = 3)

    for y_individual in fit[6]:
        a.plot(fit[-2], y_individual)

    #plot residuals
    j = 0
    ix = np.intersect1d(x, fit[-2], return_indices = True)[1]
    for i in ix:
        ysub[i] = ysub[i] - fit[0][j]
        j = j + 1
    #save peaks to file
    for i, peak in enumerate(fit[3]):
        line = str(fit[1][2 * i + 1])  + ' ' + str(np.sqrt(fit[2][2 * i + 1][2 * i + 1])) + ' ' + str(peak) + ' ' + str(fit[4][i]) + ' ' + str(fit[1][-3]) + ' ' + str(np.sqrt(fit[2][-3][-3])) + ' ' + str(fit[1][-2]) + ' ' + str(fit[1][-1]) + '\n'
        file.write(line)

file.close() 

a.plot(x, ysub)
plt.show()



