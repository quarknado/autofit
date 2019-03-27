## autosignal.py
## Signal processing functions for automatic fitting of transfer reaction spectra
## Ben Cropper 2019

import scipy as sp
import numpy as np
import pylab as pl
from scipy.signal import find_peaks_cwt
from scipy.signal import gaussian
import glob
import pandas as pd
import matplotlib.pyplot as plt

import adminfunctions as ad

smoothlength = 200 #number of points in smooting gaussian

def smoothe(yvals, smoothwidth):
    #need to smooth over the noise so I can more effectively grab fitting regions
    g = gaussian(smoothlength, smoothwidth)
    g = g/np.trapz(g) #normalise the gaussian so the scale is preserved after convolution

    return np.convolve(yvals,g,'same')

def getthreshold():
    while True:    
        try:
            thresh = float(input('input a threshold where regions above which are considered groups of peaks'))
            if thresh < 0: raise ValueError
            break
        except ValueError:
            print('invalid threshold entered')
    return thresh

def clipspectrum(xvalues, yvalues, smoothedy, threshold):
    #get where the smoothed spectrum goes above the threshold
    yix = np.where(smoothedy > threshold)

    #iplement these positions on the main spectrum
    xthreshold = xvalues[yix]
    ythreshold = yvalues[yix]
    #pl.plot(xthreshold,ythreshold, '.')
    #pl.show()

    return xthreshold, ythreshold

def split_spectrum(xthresh, ythreshold, min_region_size):

    #now separate these into different objects to be fitted separately
    peak_regions = []

    #have these temporary lists to push x and y values onto
    xreg = []
    yreg = []

    #go along the spectrum and separate it into peak regions by detecting where the gaps are 
    for i in range(len(xthresh)):
        xreg.append(xthresh[i])
        yreg.append(ythreshold[i])
    
    
        if i != len(xthresh) - 1: #checking we aren't at the end of the array
            if xthresh[i] != xthresh[i+1] - 1:
                if len(xreg) > min_region_size: #make sure these regions are actually big enough to fit, wider than the template peak width
                    peak_regions.append([xreg,yreg])
                #clean these temporary things up
                xreg = []
                yreg = []

        else:
            if len(xreg) > min_region_size:
                peak_regions.append([xreg,yreg])

    return peak_regions

def plotall(regions, xvalues, yvalues):

    fig = pl.figure(figsize = (19, 10))
    nrows = np.ceil(len(regions)/4)
    ncols = 8
    i = 1
    for x,y in regions:
        ax1 = fig.add_subplot(nrows,ncols,i)
        ax1.plot(x,y)   

        #next find where it is on the spectrum. This is so very small regions can be
        #shown what they look like in the bigger picture
        centroid = int(np.mean(x))
        x2, y2 = context(xvalues, yvalues, centroid, 100)
    
        ax2 = fig.add_subplot(nrows, ncols, i+1)
        ax2.plot(x2,y2, color = 'xkcd:orange')
        i = i+2

    fig.tight_layout()
    pl.show()

def plotcompare(xpeak, ypeak, xspectrum, yspectrum):
        fig = pl.figure()
        ax1 = fig.add_subplot(121)
        ax1.plot(xpeak, ypeak) #just plot out the fitting region

        #next find where it is on the spectrum. This is so very small regions can be
        #shown what they look like in the bigger picture
        centroid = int(np.mean(x))
        x2, y2 = context(xspectrum, yspectrum, centroid, 100)
    
        ax2 = fig.add_subplot(122)
        ax2.plot(x2,y2)
        pl.show()


def context(xval, yval, cent, window):
    cix = np.where(np.ceil(xval) == cent )[0][0]#trying to get everything as a whole number
        
    #find the indices 100 channels either side of the centroid, or less if it reaches an end to the spectrum     
    try:
        xix = np.arange(cix - window, cix + window)
        exceptor = xval[xix]
        if cix < window:
            raise IndexError
    except IndexError:
        if cix < window:
            xix = np.arange(0, cix + window)
        else:
            xix = np.arange(cix - window, len(xval))

    #now visualise this    
    #print(xix)
    xcontext = xval[xix]
    ycontext = yval[xix]

    return xcontext, ycontext
   
    
def boundconvert(x,y,ub,lb):
    lb2 = np.where(x > lb)
    ub2 = np.where(x < ub)

    inds = np.intersect1d(lb2,ub2)
    ub2 = np.max(inds)
    lb2 = np.min(inds)
    return [ub2,lb2]

def contaminant_clipper(x, yinit):    
    yclip = yinit

    while True:
        
        yclipbuffer = yclip
        
        clip_spectrum = ad.spectrum_plotter(x, yclip)
        clip_spectrum.show()        
        
        #get the user to input a valid set of bounds in which to delete the contaminant
        while True:    
            lb = input('input a lower bound on the contaminant in channels')
            ub = input('input an upper bound on the contaminant in channels')
            if np.chararray.isdigit(ub + lb):
                lb = int(lb)
                ub = int(ub)
                if (lb > min(x) and ub > lb and ub < max(x)):
                    break             
            print('Invalid bounds, try again')

        #ub, lb = boundconvert(x,y,ub,lb)
        x21 = np.where(x < lb)[0]
        x22 = np.where(x > ub)[0]    
    
        y21 = yclip[x21]
        y22 = yclip[x22]

        length_difference  = len(x) - len(x21) - len(x22)

        y20 = np.zeros(length_difference)


        ynew = np.append(y21, y20)
        ynew = np.append(ynew, y22)

        yclip = ynew
        
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


    return ynew



















