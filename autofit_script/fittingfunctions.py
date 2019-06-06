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
def lnlike(p0, x, y, background, mets = None):
    n_end = 1    
    if background:
        n_end = n_end + 2
    if mets == None:
        n_end = n_end + 2

    # get model for these parameters:
    
    #need to build up a model from multiple peaks in p0
    ymod = np.zeros(len(y))
    if mets == None: #in this case, all parameters should be fit, including the shape parameters
        #p0 goes amplitude, position, amplitude, ......, position, width, r, beta
        # so the number of peaks is the length of p0, -3 for the shape parameters, /2 for the 2 extra parameters for each peak.
        npeaks = int((len(p0)- n_end)/2)
        for i in range(npeaks):
            #add each peak in turn 
            #    AMPLITUDE   POSITION       WIDTH   R       BETA  
            p1 = [p0[i * 2], p0[i * 2 + 1], p0[-3], p0[-2], p0[-1]]
            #add the function of these parameters            
            ymod += gf3(p1,x)
    else:
        #if r and beta fixed, there's 2 fewer parameters
        npeaks = int((len(p0)-n_end)/2)
        for i in range(npeaks):
            #here r and beta are from the meta parameters
            p1 = [p0[i * 2], p0[i * 2 + 1], p0[-1], mets[1], mets[2]]
            ymod += gf3(p1,x)
    
    if background:
        ymod += p0[-n_end] * x + p0[-n_end + 1]

    # Poisson loglikelihood for the model compared to the data:
    try:
        ll = np.sum(ymod[np.where(y!=0)]*np.log(y[np.where(y!=0)])) - np.sum(y) - np.sum(ymod[np.where(ymod!=0)]*np.log(ymod[np.where(ymod!=0)]) - ymod[np.where(ymod!=0)])
    except TypeError:
        print('ymod = ', ymod)
        print('y = ', y)
               
    return ll



nll = lambda *args: -lnlike(*args) #lambda function just returning -log likelihood

def pyield(A, s, r, b):

    part1 = s * np.sqrt(2 * np.pi) * A * (1 - r/100)
    part2 = 2 * A * r * b/100 * np.exp(-(s**2)/(2 * b**2))
    
    return part1 + part2

def fit(x, y, muarr, sig, FWHM, r = 50, beta = None, rbfix = True, background = False, Aarr = None, fig = True, posfix = False):
    
    if beta == None: beta = FWHM

    if Aarr is None:
        peakix = np.intersect1d(x, muarr, return_indices = True)[1].astype(int)
        Aarr = np.take(y, peakix)
        

    #default bounds for A and mu
    #A can't be less than 0 (yay the upside down peaks that you sometimes get in gf3 are no more)
    #the peak should be on the bit of spectrum that you're fitting, so mu bounded with x
    if not posfix: bnd = ((0.,None), (min(x),max(x)), )

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
        try:        
            p = [pos, Aarr[i]]
        except IndexError:
            print('\n\n\n',muarr,'\n\n\n', Aarr)
            raise IndexError
        peakparams.append(p)
    

    p0 = []
    for pos, amp in peakparams:
        p0.append(amp)
        p0.append(pos)
    
    yieldarr = []
    yerrarr = []
    bnds = []
    for peak in muarr:
        if posfix: bnd = ((0.,None), (peak - 0.1, peak + 0.1 ))
        for bound in bnd:
            bnds.append(bound)

    nopeaks = len(muarr)
    ymod = 0

    if fig:
        fitplot = plt.figure()
        fitax = fitplot.add_subplot(111)
    individual_peaks = []

    if background:
        bg_bnd = ((None, 0.), (0., max(y)))

        for bnd in bg_bnd:
            bnds.append(bnd)
        
        for i in range(2): p0.append(0) #default 0 background

    if not rbfix:
        for bnd in metabnd:
            bnds.append(bnd)
        for met in mets:
            p0.append(met)
        result = op.minimize(nll, p0, bounds=bnds, args=(x, y, background))

    else:
        bnds.append(metabnd[0])
        p0.append(mets[0])
        result = op.minimize(nll, p0, bounds=bnds, args=(x,y, background, mets))
        
    p1 = result["x"]    
    p1cov = result["hess_inv"].todense()

    p1 = p1_fixer(p1, rbfix, background, mets)
    p1cov = cov_fixer(p1cov, rbfix, background)            

    for i in range(nopeaks):
        p2 = [p1[i * 2], p1[i * 2 + 1], p1[-3], p1[-2], p1[-1]]
        thispeak = gf3(p2,x)            
        ymod += thispeak
        if fig: fitax.plot(x,thispeak)
        individual_peaks.append(thispeak)

        yiel = pyield(p2[0],p2[2],p2[3],p2[4])
        yieldarr.append(yiel)

        yielerr = yerr(yiel,p2[0],p2[2],p2[3],p2[4], shrinkcov(p1cov, i)) 
        yerrarr.append(yielerr)

    ymod += p1[-5] * x + p1[-4]
    if fig: fitax.plot(x,p1[-5] * x + p1[-4])       


    #print('yield = ', yieldarr, ' +- ', yerrarr,)
    if fig: fitax.plot(x,y, 'b')
    if fig: fitax.plot(x,ymod, 'r', alpha = 0.7,)

    if not fig: fitplot = None

    return(ymod, p1, p1cov, yieldarr, yerrarr, fitplot, individual_peaks, x, y)
    
def yerr(Y,A,s,r,b,cov):
    
    dYdA = Y/A
    dYdr = (Y - s * np.sqrt(2 * np.pi) * A)/r
    dYds = np.sqrt(2 * np.pi) * A * (1 - r/100) - (2 * A * r * s * np.exp(-s**2/(2 * b**2)))/ (100 * b)
    dYdb = (2/100) * A * r *  np.exp(-s**2/(2 * b**2)) * (1 + s**2/b**2) 


    vAA, vAm, vAg, vao, vAs, vAr, vAb = cov[0]
    vmA, vmm, vmg, vmo, vms, vmr, vbb = cov[1]
    vgA, vgm, vgg, vgo, vgs, vgr, vgb = cov[2]   
    voA, vom, vog, voo, vos, vor, vob = cov[3]    
    vsA, vsm, vsg, vso, vss, vsr, vsb = cov[4]
    vrA, vrm, vrg, vro, vrs, vrr, vrb = cov[5]
    vbA, vbm, vbg, vbo, vbs, vbr, vbb = cov[6]
    
    #print(cov)

    vY = (dYdA**2 * vAA) + (dYds**2 * vss) + (dYdr**2 * vrr) + (dYdb**2 * vbb) + (2 * dYdA * dYds * vAs) + (2 * dYdA * dYdr * vAr) + (2 * dYdA * dYdb * vAb) + (2 * dYds * dYdr * vsr) + (2 * dYds * dYdb * vsb) + (2 * dYdr * dYdb * vrb)
    return(np.sqrt(vY))


def shrinkcov(cov,i):
    #print(cov)
    deletelist = []
    for j,row in enumerate(cov):
        if j != 2*i:
            if j != 2*i +1:
                if not np.array_equal(row, cov[-1]):
                    if not np.array_equal(row, cov[-2]):
                        if not np.array_equal(row, cov[-3]):
                            if not np.array_equal(row, cov[-4]):
                                if not np.array_equal(row, cov[-5]):
                                    deletelist.append(j)


    cov = np.delete(cov,deletelist, axis = 0)
    cov = np.delete(cov,deletelist, axis = 1)
    #print('\n', cov)
    return(cov)
        
def p1_fixer(p1, rbfix, background, mets):
    #want to 'fix' p1 so it reads amp, pos, amp, pos,....,offset, grad, sig, r, beta

    #how many parameters are currently on the 'end' of p1?
    #p1 is currently the list of parameters that have been fit
    n_end = 1    
    if background: n_end = n_end + 2
    if not rbfix: n_end = n_end + 2
    
    #by default just the width has to be there, let's get it. will be -3 if rb not fixed, -1 if they are
    if rbfix:
        sig, r, beta = [p1[-1], mets[1], mets[2]]
    else:
        sig, r, beta = [p1[-3], p1[-2], p1[-1]]

    if background:
        grad, off = [p1[-n_end], p1[-n_end + 1]]
    else:
        grad, off = [0, 0]

    p1 = p1[0:-n_end]

    p1 = np.append(p1, [grad, off, sig, r, beta])

    return(p1)

def cov_fixer(cov, rbfix, background):    
    if rbfix:
        if not background:
            scov = cov[-1]
            scovlong = np.insert(scov, -1, [0.,0.])
            scovlong = np.expand_dims(scovlong, axis = 1)
            scovshort = np.delete(scov, -1)
            scovshort = np.append(scovshort, [0, 0])
            scovshort = np.expand_dims(scovshort, axis = 0)
            cov = np.delete(cov,-1, axis = 0)
            cov = np.delete(cov,-1, axis = 1)
            cov = add2(cov)
            cov = np.append(cov, scovshort, axis = 0)
            cov = np.append(cov, scovlong, axis = 1)
 
        cov = add2(cov)
    else:
        if not background:
            scov = cov[-3]
            rcov = cov[-2]
            bcov = cov[-1]

            scovlong = np.insert(scov, -3, [0.,0.])
            rcovlong = np.insert(rcov, -3, [0.,0.])
            bcovlong = np.insert(bcov, -3, [0.,0.])

            scovlong = np.expand_dims(scovlong, axis = 1)
            rcovlong = np.expand_dims(rcovlong, axis = 1)
            bcovlong = np.expand_dims(bcovlong, axis = 1)


            scovshort = np.delete(scov, -1)
            scovshort = np.delete(scovshort, -1)
            scovshort = np.delete(scovshort, -1)
            rcovshort = np.delete(rcov, -1)
            rcovshort = np.delete(rcovshort, -1)
            rcovshort = np.delete(rcovshort, -1)
            bcovshort = np.delete(bcov, -1)
            bcovshort = np.delete(bcovshort, -1)
            bcovshort = np.delete(bcovshort, -1)

            scovshort = np.append(scovshort, [0, 0])
            rcovshort = np.append(rcovshort, [0, 0])
            bcovshort = np.append(bcovshort, [0, 0])

            scovshort = np.expand_dims(scovshort, axis = 0)            
            rcovshort = np.expand_dims(rcovshort, axis = 0)
            bcovshort = np.expand_dims(bcovshort, axis = 0)
            

            cov = np.delete(cov,-1, axis = 0)
            cov = np.delete(cov,-1, axis = 0)
            cov = np.delete(cov,-1, axis = 0)
            cov = np.delete(cov,-1, axis = 1)
            cov = np.delete(cov,-1, axis = 1)
            cov = np.delete(cov,-1, axis = 1)


            cov = add2(cov)


            cov = np.append(cov, scovshort, axis = 0)
            cov = np.append(cov, rcovshort, axis = 0)
            cov = np.append(cov, bcovshort, axis = 0)

            cov = np.append(cov, scovlong, axis = 1)
            cov = np.append(cov, rcovlong, axis = 1)
            cov = np.append(cov, bcovlong, axis = 1)
            
        else: pass

    return(cov)

def add2(cov):
    zeros = np.zeros(len(cov) * 2).reshape(2,len(cov))
    zeros2 = np.zeros((len(cov)+2) * 2).reshape(len(cov) + 2,2)
    cov2 = np.concatenate((cov, zeros), axis = 0)
    cov3 = np.concatenate((cov2, zeros2), axis = 1)
    return(cov3)
