## autosignal.py
## Signal processing functions for automatic fitting of transfer reaction spectra
## Ben Cropper 2019

import scipy as sp
import numpy as np
from scipy.signal import find_peaks_cwt
from scipy.signal import gaussian
import matplotlib.pyplot as plt
import random

import adminfunctions as ad

#number of points in smoothing gaussian. This is fairly abitrary, it's just to make sure all of the gaussian is recorded, and to a decent resolution
#it's to convolve the whole spectrum with a gaussian equal to the width of the FWHM. It should be shorter than the length of the spectrum, considerably longer than the FWHM
smoothlength = 200

def baseline_subtraction(x,y):    
    
    while True:    
        try:        
            base_split = int(input('How many pieces should the spectrum be separated into?\nThe channel with the minimum counts in each segment will be a point on a segmented linearly interpolated background.'))
            if base_split < 0: raise ValueError
            break
        except ValueError:
            print('Invalid number of sections entered')

    x_split_list = np.array_split(x,base_split)
    y_split_list = np.array_split(y,base_split)

    min_x_list = []

    for i, partx in enumerate(x_split_list):
        
        min_index = np.where(y_split_list[i] == min(y_split_list[i]))[0]        
        #sometimes there will be more than one channel with the minimum y. As to not bias this algortithm, I'll pick this randomly.        
        if len(min_index) > 1:
            min_index = random.choice(min_index)
        min_index = int(min_index)
        min_x = partx[min_index]
        min_x_list.append(min_x)

    min_indices = np.intersect1d(x, min_x_list, return_indices = True)[1]
    #print(min_indices)

    min_y_list = y[min_indices]

    interp_function = sp.interpolate.interp1d(min_x_list, min_y_list, kind = 'cubic', fill_value = 'extrapolate')    

    plt.plot(x,y)
    plt.plot(x, interp_function(x))
    plt.plot(x, y - interp_function(x))

    plt.show()

    redo = ad.y_or_n('Would you like to reselect the boundaries?')
    if redo == 'y':
        return baseline_subtraction(x,y)
    
    return( y - interp_function(x))

    

        

def smoothe(yvals, smoothwidth, length = smoothlength):
    #need to smooth over the noise so I can more effectively grab fitting regions
    #this is essentially a low-pass filter
    g = gaussian(length, smoothwidth)
    g = g/np.trapz(g) #normalise the gaussian so the scale is preserved after convolution. trapz is trapezium rule integration

    return np.convolve(yvals,g,'same')

def getthreshold():
    '''
    quality of life function to prompt the user to input a threshold. It makes the user input a positive number as a threshold 
    '''
    while True: #pins the user in a loop until they get it right    
        try:
            #errors raised if the input can't be converted to a float of if the number is less than 0
            thresh = float(input('input a threshold where regions above which are considered groups of peaks'))
            if thresh < 0: raise ValueError
            break
        except ValueError:
            print('invalid threshold entered')
    return thresh

def clipspectrum(xvalues, yvalues, smoothedy, threshold):
    '''
    returns the parts of a spectrum where the y values exceed a threshold, deleting thee rest
    '''
    #get where the smoothed spectrum goes above the threshold
    yix = np.where(smoothedy > threshold)

    #return the values at these positions on the main spectrum
    xthreshold = xvalues[yix]
    ythreshold = yvalues[yix]
    #plt.plot(xthreshold,ythreshold, '.')
    #plt.show()

    return xthreshold, ythreshold

def split_spectrum(xthresh, ythreshold, min_region_size):
    '''
    this function gets arrays that have gaps in them and splits them into separate arrays. In other words, if you have 2 peaks above the threshold and clipspectrum()
    has deleted everything in between, this function returns a list of x and y arrays separately i.e. [[xarr1, yarr1], [xarr2, yarr2]....] 
    '''
    #define the master list of peak regions
    peak_regions = []

    #have these temporary lists to push x and y values onto
    xreg = []
    yreg = []

    #go along the spectrum and separate it into peak regions by detecting where the gaps are 
    for i in range(len(xthresh)):
        #firstly, stick the ith element of x and y on their list
        xreg.append(xthresh[i])
        yreg.append(ythreshold[i])
    
    
        if i != len(xthresh) - 1: #checking we aren't at the end of the array
            if xthresh[i] != xthresh[i+1] - 1: #if there's a gap in x, this is where we want the peak region to end
                if len(xreg) > min_region_size: #make sure these regions are actually big enough to fit, wider than the template peak width
                    xreg = np.array(xreg)
                    yreg = np.array(yreg)
                    peak_regions.append([xreg,yreg])
                #clean these temporary things up
                xreg = []
                yreg = []

        else: #this is for if we've reached the end of the spectrum - save the region currently being created and the loop ends
            if len(xreg) > min_region_size:
                xreg = np.array(xreg)
                yreg = np.array(yreg)
                peak_regions.append([xreg,yreg])

    return peak_regions

def plotall(regions, xvalues, yvalues):
    '''
    this one is to plot all of the peak regions and their contexts on one figure, so the user can evaluate whether these are correct
    '''
    fig = plt.figure(figsize = (19, 10))
    #set it up so that there are 8 columns (containing 4 regions, 2 plots each)
    nrows = np.ceil(len(regions)/4)
    ncols = 8
    i = 1 #can't just use enumerate for regions because index skips two each time
    for x,y in regions:
        #plot the fitting region        
        ax1 = fig.add_subplot(nrows,ncols,i)
        ax1.plot(x,y)   

        #next find where it is on the spectrum. This is so very small regions can be
        #shown what they look like in the bigger picture
        centre = int(np.mean(x))
        x2, y2 = context(xvalues, yvalues, centre, 100)
        
        #now plot the fitting region context
        ax2 = fig.add_subplot(nrows, ncols, i+1)
        ax2.plot(x2,y2, color = 'xkcd:orange')
        i = i+2

    fig.tight_layout()
    fig.show()

def context(xval, yval, cent, window):
    '''
    Plots the context of a fitting region. This is useful when you don't see much from the plot of the fitting region itself.
    Plotting 100 channels either side helps you know what's going on with a potentially spurious peak.
    '''
    cix = np.where(np.ceil(xval) == cent )[0][0]#trying to get the centroid indices as a whole number
        
    #find the indices 100 channels either side of the centre, or less if it reaches an end to the spectrum     
    try:
        xix = np.arange(cix - window, cix + window)
        exceptor = xval[xix]#this will raise an index error if the centre is less than 'window' channls away from the edge
        if cix < window: #raises the exception if the middle is close to the '0' end of the spectrum
            raise IndexError
    except IndexError: #truncates the window at the end or beginning of the spectrum
        if cix < window:
            xix = np.arange(0, cix + window)
        else:
            xix = np.arange(cix - window, len(xval))

    
    xcontext = xval[xix]
    ycontext = yval[xix]

    return xcontext, ycontext
   

def contaminant_clipper(x, yinit, ybase):

    '''
    clips bits out of the spectrum after prompting the user to do so
    '''    
    yclipfit = yinit #clip things out of this yclip variable so yinit doesn't change
    yclipfind = ybase

    while True:
        
        #set a buffer such that it can be reinstated to the last one if the user changes their mind
        yclipfitbuffer = yclipfit
        yclipfindbuffer = yclipfind                

        #plot what the current spectrum is
        clip_spectrum = ad.spectrum_plotter(x, yclipfind)
        clip_spectrum.show()        
        
        #get the user to input a valid set of bounds in which to delete the contaminant
        while True:    
            lb = input('input a lower bound on the contaminant in channels')
            ub = input('input an upper bound on the contaminant in channels')
            if np.chararray.isdigit(ub + lb): #has to all be digits
                lb = int(lb)
                ub = int(ub)
                if (lb > min(x) and ub > lb and ub < max(x)): #lower and upper bounds have to be in the spectrum, upper bound has to be bigger than lower
                    break             
            print('Invalid bounds, try again')

        #get which indices x is less than the lower bound and greater than the upper bound set
        x21 = np.where(x < lb)[0]
        x22 = np.where(x > ub)[0]    
        
        #get y arrays above and below the contaminant
        y21fit = yclipfit[x21]
        y22fit = yclipfit[x22]

        y21find = yclipfind[x21]
        y22find = yclipfind[x22]

        #how many channels does the contaminant have
        length_difference  = len(x) - len(x21) - len(x22)

        #get an array of 0s to replace the contaminant ys with
        y20 = np.zeros(length_difference)

        #now to stick together the arrays
        #need the ys below the contaminant, then the zeroes where the contaminant was, and then the ys above the contaminant
        ynewfit = np.append(y21fit, y20)
        ynewfit = np.append(ynewfit, y22fit)

        ynewfind = np.append(y21find, y20)
        ynewfind = np.append(ynewfind, y22find)

        #now replace the old ys with the new one
        yclipfit = ynewfit
        yclipfind = ynewfind
        
        #now plot that and ask if they're happy with it
        clip_spectrum = ad.spectrum_plotter(x, yclipfind)
        clip_spectrum.show() 
        
        confirmclip = ad.y_or_n('Are you happy with this removal?')

        #what to do now? if they aren't happy, ask if they want another go. if they think there aren't any contaminants, return to yinit
        #if there still are, then they can have another go
        if confirmclip == 'n':

            cancel = ad.y_or_n('Are you sure there are any contaminants?')
            if cancel == 'n':
                yclipfit = yinit
                yclipfind = ybase
                break
            if cancel == 'y':
                yclipfit = yclipfitbuffer
                yclipfind = yclipfindbuffer
                continue
            
        else: #all good, so can continue
            pass

        #do they want to get rid of more contaminants? continue loop if they do, don't if they don't
        moreclip = ad.y_or_n('Would you like to remove any more contaminants?')
        if moreclip == 'y':
            continue
        else:
            break


    return yclipfit, yclipfind

def peak_finder(peak_regions, w, fwhm):
    '''
    Tries to find the peak positions for each region, first by applying a low-pass filter
    then  a continuous wavelet transform, varying wavelet width over the plausible peak widths and finding ridges of local maxima in the cwt matrix.
    (see the scipy docs for a better description of how this works. It is very interesting!)
    '''
    regions_positions = []
    peak_positions = []
    for x2, y in peak_regions:
        #low pass filter. Narrow gaussian, so only smoothes a little bit.
        ys = smoothe(y, w/4, length = len(y))
    
        #now do the cwt. the widths here are the widths of the wavelet you convolve the signal width
        #A peak which persists when you convolve the signal with wavelets of different widths will be considered a peak
        reg_positions = find_peaks_cwt(ys, widths = np.arange( w ,fwhm))
        #this one finds the absolute peak positions on the spectrum and lists them.
        peak_positions += (min(x2) + find_peaks_cwt(ys, widths = np.arange( w ,fwhm))).tolist()
        #this one sorts them into regions
        regions_positions.append(min(x2) + reg_positions)
    return peak_positions, regions_positions

















