## fittingfunctions.py
## Fitting functions for automatic fitting of transfer reaction spectra
## Ben Cropper 2019

import numpy as np
from scipy.signal import find_peaks_cwt
from scipy.special import erfc
import matplotlib.pyplot as plt
import scipy.optimize as op

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
def lnlike(p0, x, y, mets = None):
    # get model for these parameters:

    #need to build up a model from multiple peaks in p0
    ymod = 0
    if mets == None: #in this case, all parameters should be fit, including the shape parameters
        #p0 goes amplitude, position, amplitude, ......, position, width, r, beta
        # so the number of peaks is the length of p0, -3 for the shape parameters, /2 for the 2 extra parameters for each peak.
        npeaks = int((len(p0)-3)/2)
        for i in range(npeaks):
            #add each peak in turn 
            #    AMPLITUDE   POSITION       WIDTH   R       BETA  
            p1 = [p0[i * 2], p0[i * 2 + 1], p0[-3], p0[-2], p0[-1]]
            #add the function of these parameters            
            ymod += gf3(p1,x)
    else:
        #if r and beta fixed, there's 2 fewer parameters
        npeaks = int((len(p0)-1)/2)
        for i in range(npeaks):
            #here r and beta are from the meta parameters
            p1 = [p0[i * 2], p0[i * 2 + 1], p0[-1], mets[1], mets[2]]
            ymod += gf3(p1,x)

    # Poisson loglikelihood for the model compared to the data:
    ll = np.sum(ymod[np.where(y!=0)]*np.log(y[np.where(y!=0)])) - np.sum(y) - np.sum(ymod[np.where(ymod!=0)]*np.log(ymod[np.where(ymod!=0)]) - ymod[np.where(ymod!=0)])
    return ll



nll = lambda *args: -lnlike(*args) #lambda function just returning -log likelihood

def pyield(A, s, r, b):

    part1 = s * np.sqrt(2 * np.pi) * A * (1 - r/100)
    part2 = 2 * A * r * b/100 * np.exp(-(s**2)/(2 * b**2))
    
    return part1 + part2

def fit(x, y, muarr, sig, FWHM, r = 50, beta = None, rbfix = True):
    
    if beta == None: beta = FWHM
    
    peakix = np.intersect1d(x, muarr, return_indices = True)[1].astype(int)
    Aarr = np.take(y, peakix)
        


    #default bounds for A and mu
    #A can't be less than 0 (yay the upside down peaks that you sometimes get in gf3 are no more)
    #the peak should be on the bit of spectrum that you're fitting, so mu bounded with x
    bnd = ((0.,None), (min(x),max(x)), )
    #bounds for metaparameters sigma, r, and beta (they are 'meta' because they should remain constant across peaks locally)
    #set sigma bound at the resolution of the detector - it's about 4 channels for Munich using the binning that I am
    #This can change though. I've set the upper bound to the FWHM input earlier just so it doesn't try to fit the background
    #It can't be negative either
    #R can only be between 0 and 100 since it's a percentage of the height
    #Beta can be any positive number
    metabnd = ((0.1,FWHM), (0.1, 100), (0.1,None))
    mets = [sig, r, beta]


    peakparams = []
    for i, pos in enumerate(muarr):
        p = [pos, Aarr[i]]
        peakparams.append(p)
    

    p0 = []
    for pos, amp in peakparams:
        p0.append(amp)
        p0.append(pos)
    
    yieldarr = []
    yerrarr = []
    bnds = []
    for peak in muarr:
        for bound in bnd:
            bnds.append(bound)

    nopeaks = len(muarr)
    ymod = 0

    fitplot = plt.figure()
    fitax = fitplot.add_subplot(111)
    individual_peaks = []

    if rbfix == False:
        for bnd in metabnd:
            bnds.append(bnd)
        for met in mets:
            p0.append(met)
        result = op.minimize(nll, p0, bounds=bnds, args=(x, y))
        p1 = result["x"]
        p1cov = result["hess_inv"].todense()
        
        for i in range(nopeaks):
            p2 = [p1[i * 2], p1[i * 2 + 1], p1[-3], p1[-2], p1[-1]]
            thispeak = gf3(p2,x)            
            ymod += thispeak
            fitax.plot(x,thispeak)
            individual_peaks.append(thispeak)
       
            yiel = pyield(p2[0],p2[2],p2[3],p2[4])
            yieldarr.append(yiel)

            yielerr = yerr(yiel,p2[0],p2[2],p2[3],p2[4], shrinkcov(p1cov, i, rbfix)) 
            yerrarr.append(yielerr)

    else:
        bnds.append(metabnd[0])
        p0.append(mets[0])
        result = op.minimize(nll, p0, bounds=bnds, args=(x,y, mets))
        p1 = result["x"]
        p1cov = result["hess_inv"].todense()
        p1 = np.concatenate((p1, np.array([mets[1], mets[2]])), axis = None) 
        p1cov = zerorb(p1cov)

       

        for i in range(nopeaks):
            p2 = [p1[i * 2], p1[i * 2 + 1], p1[-3], mets[1], mets[2]]
            thispeak = gf3(p2,x)            
            ymod += thispeak
            fitax.plot(x,thispeak)
            individual_peaks.append(thispeak)

            yiel = pyield(p2[0],p2[2],p2[3],p2[4])
            yieldarr.append(yiel)

            yielerr = yerr(yiel,p2[0],p2[2],p2[3],p2[4], shrinkcov(p1cov, i, rbfix)) 
            yerrarr.append(yielerr)

       


    #print('yield = ', yieldarr, ' +- ', yerrarr,)
    fitax.plot(x,y, 'b')
    fitax.plot(x,ymod, 'r', alpha = 0.7,)

    return(ymod, p1, p1cov, yieldarr, yerrarr, fitplot, individual_peaks, x, y)
    
def yerr(Y,A,s,r,b,cov):
    
    dYdA = Y/A
    dYdr = (Y - s * np.sqrt(2 * np.pi) * A)/r
    dYds = np.sqrt(2 * np.pi) * A * (1 - r/100) - (2 * A * r * s * np.exp(-s**2/(2 * b**2)))/ (100 * b)
    dYdb = (2/100) * A * r *  np.exp(-s**2/(2 * b**2)) * (1 + s**2/b**2)
    
    vAA, vAm, vAs,vAr, vAb = cov[0]
    vmA, vmm, vms,vmr, vbb = cov[1]
    vsA, vsm, vss,vsr, vsb = cov[2]
    vrA, vrm, vrs,vrr, vrb = cov[3]
    vbA, vbm, vbs,vbr, vbb = cov[4]
    
    #print(cov)

    vY = (dYdA**2 * vAA) + (dYds**2 * vss) + (dYdr**2 * vrr) + (dYdb**2 * vbb) + (2 * dYdA * dYds * vAs) + (2 * dYdA * dYdr * vAr) + (2 * dYdA * dYdb * vAb) + (2 * dYds * dYdr * vsr) + (2 * dYds * dYdb * vsb) + (2 * dYdr * dYdb * vrb)
    return(np.sqrt(vY))


def shrinkcov(cov,i, fixrb):
    #print(cov)
    deletelist = []
    for j,row in enumerate(cov):
        if j != 2*i:
            if j != 2*i +1:
                if not np.array_equal(row, cov[-1]):
                    if not np.array_equal(row, cov[-2]):
                        if not np.array_equal(row, cov[-3]):
                            deletelist.append(j)


    cov = np.delete(cov,deletelist, axis = 0)
    cov = np.delete(cov,deletelist, axis = 1)
    #print('\n', cov)
    return(cov)
        

def zerorb(cov):
    zeros = np.zeros(len(cov) * 2).reshape(2,len(cov))
    zeros2 = np.zeros((len(cov)+2) * 2).reshape(len(cov) + 2,2)
    cov2 = np.concatenate((cov, zeros), axis = 0)
    cov3 = np.concatenate((cov2, zeros2), axis = 1)
    return(cov3)
