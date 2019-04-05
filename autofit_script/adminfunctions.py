## adminfunctions.py
## Functions for admin purposes (reading files, getting yes or no answers etc) for automatic fitting of transfer reaction spectra
## Ben Cropper 2019

import os
import numpy as np
import matplotlib.pyplot as plt

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

def y_or_n(question):
    while True:        
        yorn = input(question + 'y/n')
        if yorn == 'y': break
        elif yorn == 'n' : break 
        else: print('invalid input, please try again')
        
    return yorn

def listfiles(directory):
    #set an index, since it's going to print out the index before each file so the user can
    #select which file to access more easily
    fileindex = 0
    #simple loop through the directory
    for roots, dirs, filenames in os.walk(directory):
        for f in filenames:
            fileindex = fileindex + 1
            #need to print the name of each file out so the user can select which one to open 
            print( str(fileindex) + ': ' + f)
    #returns the total number of files so the program can detect whether the
    #selected file number is valid (ie falls in between 0 and the number of files)
    return fileindex

def inputfileselector(fileno): #fileno is the number of files in the directory
    #Initialize the desired quantity
    fileselect = None
    while fileselect == None:
        #prompt the user to input their selection
        fileselect = input('Which file would you like to open?\n')
        try:
            
            #failure test
            #in this case, it checks to see if the input matches a file index
            
            #create a list of the possible file indices
            fileindices = list(range(1, fileno + 1))
            #filetest doesn't do annything, but if there's an error in the assignment, it
            #triggers the error
            filetest = fileindices.index(int(fileselect))

        except:
            fileselect = None
    return int(fileselect)

def openfilefromlist(fileselect, directory):
    for roots, dirs, filenames in os.walk(directory):
        f = filenames[fileselect-1]
    return f

def spectrum_plotter(x, y):
    fig = plt.figure(figsize = (15,7))
    ax = fig.add_subplot(111)
    ax.plot(x, y)
    ax.set_xticks(np.arange(0,max(x), 100))

    return fig

def printfit(fit, x, y):
    ymod, p1, p1cov, yieldarr, yerrarr = fit[0], fit[1], fit[2], fit[3], fit[4]
    sig, r, beta = p1[-3], p1[-2], p1[-1]
        
    nopeaks = int((len(p1) - 3)/2)
    chi2 = np.sum((y - ymod)**2/y)
    if p1cov[-1][-1] == 0:
        ndof = len(y) - len(p1) + 2
    else:
        ndof = len(y) - len(p1)
    
    print('Fit complete. Fitted channels ', min(x), ' to ', max(x), ', with ', nopeaks, ' peak(s).') 
    print('Width = ', paramprint(sig, np.sqrt(p1cov[-3][-3])), '\nR = ', paramprint(r, np.sqrt(p1cov[-2][-2])), '\nBeta = ', paramprint(beta, np.sqrt(p1cov[-1][-1])), '\nReduced chi-squared = ', chi2/ndof) 


    print('Positions and yields:')
    
    for i in range(nopeaks):
        Y = yieldarr[i]
        sY = yerrarr[i]
        pos = p1[2*i + 1]
        spos = np.sqrt(p1cov[2*i + 1][2*i + 1])
        print('[' + str(i+1) + '] ' + paramprint(pos, spos) + ', ' + paramprint(Y, sY))
        
    
def paramprint(param, error):

    try:
        errormag =  int(np.floor(np.log10(error)))
    except OverflowError: #for errors of zero    
        errormag = 0

    if errormag > -1:
        param = int(np.around(param))
        error = int(np.around(error))
    else:
        param = np.around(param, -errormag)
        error = np.around(error, -errormag)

    
    printout = str(param) + ' +- ' + str(error)
    return(printout)


