import scipy as sp
import numpy as np
import pylab as pl
from scipy.signal import find_peaks_cwt
from scipy.signal import gaussian
import glob
import pandas as pd

#Hard coded global variables
smoothlength = 200 #length of Gausian which you smoothe the spectrum with
smoothwidth = 10 #width of Gaussian which you smoothe the spectrum with
min_region_size = 5 #minimum width of the region above threshold to be considered



def file_reader(handle):
    
    print('reading ', handle)
    
    infile = open(handle,'r')
    xvals=[];ymeas=[]
    
    #need this conditional because the asciis are comma separated and txts are space separated 
    if handle[-1] == 'c':
        #this bit reads in the files
        while True:
            line = infile.readline()
            if not line: break
        
            items = line.split()
            xvals.append(float(items[0][:-1]))
            ymeas.append(float(items[1]))

    else:
        while True:
            line = infile.readline()
            if not line: break
        
            items = line.split()
            xvals.append(float(items[0]))
            ymeas.append(float(items[1]))

    infile.close()

    #x and y should be numpy arrays so we can do more stuff with them
    ymeas = np.array(ymeas)
    xvals = np.array(xvals)

    #gets rid of the 0 channel dump
    ymeas[0] = 0

    print('successfully read ', handle)

    return xvals, ymeas

def smoothe(yvals):
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

def split_spectrum(xthreshold, ythreshold):

    #now separate these into different objects to be fitted separately
    peak_regions = []

    #have these temporary lists to push x and y values onto
    xreg = []
    yreg = []

    #go along the spectrum and separate it into peak regions by detecting where the gaps are 
    for i in range(len(xthresh)):
        xreg.append(xthreshold[i])
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

def plotall(peak_regions, xvalues, yvalues):

    fig = pl.figure(figsize = (19, 10))
    nrows = np.ceil(len(peak_regions)/4)
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
        ax2.plot(x2,y2)
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
            xix = np.arange(cix - window, len(xvals))

    #now visualise this    
    #print(xix)
    xcontext = xval[xix]
    ycontext = yval[xix]

    return xcontext, ycontext

def y_or_n(question):
    while True:        
        yorn = input(question + 'y/n')
        if yorn == 'y': break
        elif yorn == 'n' : break 
        else: print('invalid input, please try again')
        
    return yorn

'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

allx = []
ally = []
class_list = []

#make list of files, which are .txt or .asc
filenames = glob.glob('*.txt')
filenames += glob.glob('*.asc')
filenames += glob.glob('*.dat')
print('The filenames are:')
print(filenames)


#loop through all of the files to get the peak regions for all of them
for f in filenames:

    xvals, ymeas = file_reader(f)

    smooth = smoothe(ymeas)
    
    #show spectrum and smoothed spectrum on the same plot
    pl.plot(xvals, ymeas)
    pl.plot(xvals, smooth)
    pl.show()


    #now to extract the peak positions
    #this needs to be on a loop to check if you're happy with the threshold
    while True:    
        #first set the threshold
        thresh = getthreshold()
        xthresh, ythresh = clipspectrum(xvals, ymeas, smooth, thresh)
        regions = split_spectrum(xthresh, ythresh)
        plotall(regions, xvals, ymeas)

        is_threshold_good = y_or_n('Are you happy with this threshold?')
        
        if is_threshold_good == 'y': break
        else: print('ok, try again with a different threshold')    


    #go through the fitting regions, display them next to their place in the wider spectrum, and prompt the user to tell if they're contaminants.
    for x,y in regions:
        plotcompare(x,y,xvals,ymeas)
        
        ispeaks = y_or_n('Does this image correspond to a good peak or set of peaks?')

        allx.append(x)
        ally.append(y)
        class_list.append(ispeaks)
        
allx = pd.Series(allx)
ally = pd.Series(ally)
class_list = pd.Series(class_list)

data = pd.concat([allx, ally, class_list], axis=1)
data.columns = ['xvalues','ymeasured','ispeaks']

data.to_csv('train_data', sep = ' ')
