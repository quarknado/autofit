## autofit.py
## Automatic fitting routine for transfer reaction spectra
## Ben Cropper 2019

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
from scipy.signal import find_peaks_cwt, find_peaks
import glob
import os
import copy
import pickle

import adminfunctions as ad
import fittingfunctions as fit
import autosignals as sig

fstgo = True

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

#yinit = sig.baseline_subtraction(x,yinit, fstgo)


'''
###########################################Contaminant Deletion##################################################
'''

cont_del =  ad.y_or_n("Would you like to remove any contaminants?")

if cont_del == 'y':
    xlist, ylist, yclip  = sig.contaminant_clipper(x,yinit,yinit, split = True)
else:
    xlist, ylist = [[x], [yinit]]
    yclip = yinit

for x2, y2 in zip(xlist, ylist):
    plt.plot(x2, y2)
plt.show()

 

'''
##########################################Threshold set and Peak Region Detection###########################################################
'''
rlist, betalist = [],[]

fitlistlist = []
totalnfits = 25
for i in range(totalnfits):  
    while True:
        if fstgo:
            spe2 = ad.spectrum_plotter(x, yclip, 5)
            spe2.show()
    
            FWHM = float(input('What is the approximate FWHM of your peaks in channels'))
            width = FWHM/(2 * np.sqrt(2*np.log(2)))
    
            len_regions = float(input('Approximately how many channels long would you like your fitting segments to be?'))
            no_regions = len(x)/len_regions     
    
        peak_regions = sig.min_region_finder(xlist, ylist, no_regions, fstgo)
        i = 0
        none_indices = []
        for x2, y2 in peak_regions:
            if len(x2) == 0:
                none_indices.append(i)         
            if fstgo: plt.plot(x2, y2)
            i = i + 1
        if fstgo: plt.show()
        peak_regions = np.delete(peak_regions, none_indices, axis = 0)
    
    
        if fstgo:
            recontaminant = ad.y_or_n('Would you like to re-clip for any more contaminants?')
        
            if recontaminant == 'y':
                yclipfit, yclipfind  = sig.contaminant_clipper(x,yclipfit, yclipfind)
                continue
            else:
                pass
    
            final_regions = ad.y_or_n('Would you like to re-tune the width and number of peaks?')
    
            if final_regions == 'y':
                continue
            else:
                break
        else: break
    




    '''
    ###########################################Fitting##################################
    '''
    if fstgo:
        bg = ad.y_or_n('Would you like to have a local linear background added to these fits?')

        if bg == 'y':
            bg = True
        else:
            bg = False

    peak_positions, region_positions = sig.peak_finder(peak_regions, width, FWHM)

    pos = np.zeros(len(x))
    pos[np.intersect1d(x, peak_positions, return_indices = True)[1]] = 1.2 * max(yclip)

    if fstgo:
        fig = plt.figure(figsize = (15,7))
        axis = fig.add_subplot(111)
        axis.set_xticks(np.arange(0,max(x), 100))
    
        axis.plot(x,yinit, color = 'xkcd:light grey')
        for x2,y in peak_regions:
            axis.plot(x2,y)
    
        axis.plot(x,pos, linewidth = 1)

        plt.show()

        template, template_pos, peak_pos = ad.template_finder(peak_regions, peak_positions, region_positions, fstgo, peaks_figure = fig)
    else:
        template, template_pos, peak_pos = ad.template_finder(peak_regions, peak_positions, region_positions, fstgo, pos = peak_pos)



    template_fit = fit.fit(template[0], template[1], template_pos, width, FWHM, rbfix = False, background = bg, fig = False)
    ad.printfit(template_fit, template[0], template[1])

    rfix, betafix = template_fit[1][-2],template_fit[1][-1]
    rlist.append(rfix)
    betalist.append(betafix)

    fitlist = []
    for i, region in enumerate(peak_regions):
        xreg = region[0]
        yreg = region[1]
        if len(region_positions[i]) == 0: continue     
        ft = fit.fit(xreg, yreg, region_positions[i], width, FWHM, r = rfix, beta = betafix, background = bg, fig = True)
        ad.printfit(ft, xreg, yreg)

        fitlist.append(ft)

    '''
    ############################################2nd pass#############################
    '''
    muarrnew = []
    aarrnew = []
    while True:
        oldlen = len(aarrnew)    
        muarrnew = []
        aarrnew = []

        for ft in fitlist:
            for i, yiel in enumerate(ft[3]):
                if not np.isnan(ft[4][i]):
                    if yiel > 0.0000001 * ft[4][i]:
                        muarrnew.append(ft[1][2 * i + 1])
                        aarrnew.append(ft[1][2 * i])
                    

    
        newlen = len(aarrnew)
    
        if oldlen == newlen: break
        fitlist = []    
        muarrnew = np.array(muarrnew)
        aarrnew = np.array(aarrnew)
    
        for i, region in enumerate(peak_regions):
            xreg = region[0]
            yreg = region[1]
            


            reg_pos = muarrnew[(muarrnew > min(xreg)) & (muarrnew < max(xreg))]
            a_arr = aarrnew[(muarrnew > min(xreg)) & (muarrnew < max(xreg))]
    

            if len(region_positions[i]) == 0: continue     
            ft = fit.fit(xreg, yreg, reg_pos, width, FWHM, r = rfix, beta = betafix, background = bg, Aarr = a_arr, fig = True)
            ad.printfit(ft, xreg, yreg)

            fitlist.append(ft)
     
    fitlistlist.append(fitlist)    
    fstgo = False

#now want to plot pos vs yield for every peak in this list.
#print(fitlistlist)

afin = []
mfin = []
smfin = []
ctrfin = []
sctrfin = []
nlist = []

for j, fitlist in enumerate(fitlistlist):
    for ft in fitlist:
        print(len(ft[3]),len(ft[6]))
        for i, yiel in enumerate(ft[3]):
            mfin.append(ft[1][2 * i + 1])
            afin.append(ft[1][2 * i])
            smfin.append(np.sqrt(ft[2][2*i + 1][2*i + 1]))
            nlist.append(j)
        for i, peak in enumerate(ft[6]):
            x2 = ft[7]
            centroid = np.sum(peak * x2)/np.sum(peak)
            ctrfin.append(centroid)
            sctrfin.append(np.sqrt(ft[2][2*i + 1][2*i + 1]))
            #sctrfin.append((np.sum(peak * x2 ** 2)/np.sum(peak) - centroid ** 2)/ np.sqrt(np.sum(peak)))

plt.errorbar(nlist, mfin, smfin, linestyle = "")
plt.show()

plt.errorbar(mfin,ctrfin,smfin, fmt = 'x')
plt.show()

poslist = []
sposlist = []
alist = []
centroidlist = []
scentroidlist = []

#print(ctrfin, '\n\n\n\n\n')



for mu, smu, a, ctr, sctr in zip(mfin, smfin, afin, ctrfin, sctrfin):
    z = 10000
    for pos, spos, h, cent, scent in zip(poslist,sposlist, alist, centroidlist, scentroidlist):
        '''p,sp = np.average(pos, weights = 1/np.array(spos)**2, returned = True)
        sp = sp ** (-0.5)
        s = np.sqrt(sp**2 + smu**2)
        z = np.abs(mu - p)/s
        '''
        c,sc = np.average(cent, weights = 1/np.array(scent)**2, returned = True)
        sc = sc ** (-0.5)
        s = np.sqrt(sc**2 + sctr**2)
        z = np.abs(ctr - c)/s
        
        if np.isnan(c): z = 10000
 
        if z < 3.5:
            cent.append(ctr)
            scent.append(sctr)
            pos.append(mu)
            spos.append(smu)
            h.append(a)
            break
    if z < 3.5: continue
    
    poslist.append([mu])
    sposlist.append([smu])
    alist.append([a])
    centroidlist.append([ctr])
    scentroidlist.append([sctr])

for i, pos in enumerate(poslist):
    plt.errorbar(np.full(len(pos), i), pos, sposlist[i], linestyle = "" )

plt.show()



muarrnew = []
aarrnew = []
#print(alist)
for pos, a in zip(poslist, alist):
    if len(pos) > 0.85 * totalnfits:
        muarrnew.append(np.average(pos))
        aarrnew.append(np.average(a))

muarrnew = np.array(muarrnew)
aarrnew = np.array(aarrnew) 

#print(muarrnew, aarrnew)
peak_regions = sig.min_region_finder(xlist, ylist, no_regions, fstgo)
i = 0
none_indices = []


for x2, y2 in peak_regions:
    if len(x2) == 0:
        none_indices.append(i)         
        i = i + 1
    peak_regions = np.delete(peak_regions, none_indices, axis = 0)



reg_pos = []
a_arr = []

for i, region in enumerate(peak_regions):
    xreg = np.array(region[0])
    yreg = np.array(region[1])
       
    reg_pos.append(muarrnew[(muarrnew > min(xreg)) & (muarrnew < max(xreg))])
    a_arr.append(aarrnew[(muarrnew > min(xreg)) & (muarrnew < max(xreg))])


  

rfix, betafix = np.average(rlist), np.average(betalist)


fitlist = []



for i, region in enumerate(peak_regions):
    xreg = np.array(region[0])
    yreg = np.array(region[1])

    if len(reg_pos[i]) == 0: continue     
    ft = fit.fit(xreg, yreg, reg_pos[i], width, FWHM, r = rfix, beta = betafix, background = bg, Aarr = a_arr[i], posfix = True)
    #ft = fit.fit(xreg, yreg, reg_pos[i], width, FWHM, rbfix = False, background = bg, Aarr = a_arr[i], posfix = True)
    ad.printfit(ft, xreg, yreg)



    fitlist.append(ft)



'''
##########################################Plot Fits and Save######################
'''

ysub = copy.deepcopy(yclip)

fitplot = plt.figure(figsize = (15,7))
a = fitplot.add_subplot(111)
a.plot(x, yinit, color = 'xkcd:grey')

ofilename = str(f[:-4]) + '_fit.txt'
fil = open(ofilename, 'w+')
topline = 'POSITION sPOSITION AREA sAREA WIDTH sWIDTH R BETA\n'
fil.write(topline)

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
        fil.write(line)

fil.close() 


pos = np.zeros(len(x))
pos_ix = []
for mu in muarrnew:
    pos_ix.append(ad.find_nearest(mu, x))
pos[pos_ix] = max(yclip) * 1.2

a.plot(x, pos, linewidth = 1)

a.plot(x, ysub)
f2 = str(f[:-4]) + '_fig.pkl'
print(os.getcwd(), f2)
pickle.dump(fitplot, open(f2, 'wb'))
plt.show()


