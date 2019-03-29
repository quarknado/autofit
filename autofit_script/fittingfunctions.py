## fittingfunctions.py
## Fitting functions for automatic fitting of transfer reaction spectra
## Ben Cropper 2019

import scipy as sp
import numpy as np
import pylab as pl
from scipy.signal import find_peaks_cwt
from scipy.signal import gaussian
import glob
import pandas as pd
import matplotlib.pyplot as plt

import autosignals as sig
import adminfunctions as ad

#Define gf3 function
def gf3(p0, x):
    #this is the fuctional form that we have been using in our manual fit
    #it is the sum of a gaussian and a skewed gaussian with the same mean
    #the extra parameters r and beta are introduced
    #r is the fraction of the height of the skewed gaussian given as a pecentage
    #beta is the 'skewneess' of the second skew gaussian
    #it is the decay constant of an exponential tail on the skewed gaussian
    #this exponential tail is convolved with a gaussian resolution function
    
    amp, mu, sigma, r, beta = p0
    
    
    #gaussian part
    ygaus = amp * (1 - r/100) * np.exp((-(x - mu)**2)/(2 * sigma**2))
    
    #'skew' gaussian part. erfc is 1 - the error function
    yskew = amp * (r/100) * np.exp((x-mu)/beta) * erfc( (x-mu)/(sigma * np.sqrt(2))  + sigma/(beta*np.sqrt(2)))
    #yskew = 0
    #ygaus = 0
    ymod = yskew + ygaus
    
    return ymod

#define log likelihood, the thing we want to maximise
def lnlike(p0, x, y):
    
    # get model for these parameters:
    npeaks = int((len(p0)-3)/2)
    ymod = 0
    for i in range(npeaks):
        p1 = [p0[i * 2], p0[i * 2 + 1], p0[-3], p0[-2], p0[-1]]
        ymod += gf3(p1,x)
    
    # Poisson loglikelihood:
    ll = np.sum(ymod[np.where(y!=0)]*np.log(y[np.where(y!=0)])) - np.sum(y) - np.sum(ymod[np.where(ymod!=0)]*np.log(ymod[np.where(ymod!=0)]) - ymod[np.where(ymod!=0)])
    
    return ll

def template_finder(peak_regions, peak_positions, reg_pos, peaks_figure):

    for i,p in enumerate(peak_positions):
        print( '[' + str(i) + ']: ' + str(p))

    template_select = None

    while template_select == None:
    
        try:
            peaks_figure.show()
            template_select = int(input('Which peak is the most suitable as a template to fix parameters to The region containing it will be what is fitted.'))
            i = 0
            for x2, y in peak_regions:
                if len(np.intersect1d(peak_positions[template_select], x2)) == 1:
                    template = [x2, y]
                    template_pos = reg_pos[i]
                i = i + 1
        except:
            template_select = None

    plt.plot(template[0], template[1])
    plt.show()

    return template, template_pos

def peak_finder(peak_regions, w, fwhm):
    regions_positions = []
    peak_positions = []
    for x2, y in peak_regions:
        ys = sig.smoothe(y, w/4, length = len(y))
    
        reg_positions = find_peaks_cwt(ys, widths = np.arange( w ,fwhm))
        peak_positions += (min(x2) + find_peaks_cwt(ys, widths = np.arange( w ,fwhm))).tolist()
        regions_positions.append(min(x2) + reg_positions)
    print(peak_positions)
    return peak_positions, regions_positions


