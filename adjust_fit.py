"""
adjust_fit.py
This is a script to be used following the execution of fitting.py.
Created on Mon 24 June 2019
@author: O. H. Wilkins
"""

import os
import pylab as plb
import matplotlib.pyplot as plt
import numpy as np
import time
import matplotlib.patches as patches
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy.optimize import least_squares
from scipy.linalg import svd
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import aplpy
import lmfit
from lmfit import minimize, Parameters, Model, fit_report, report_ci
from uncertainties import ufloat
from datetime import date,datetime

# Not all of these packages are necessary for this script, 
# but any extraneous packages in this list are needed for 
# the fitting.py script so you would have needed them anyway.

##########
# INPUTS #
##########

print 'Declaring input files.'
print ''
ParentDir = ''
ParentDate = ''       # YYYYMMDD
ParentTime = ''       # HHMM
ParentFile = ParentDir+'ALMA_fitting'+ParentDate+'_'+ParentTime

ParentTemp = ParentFile+'_Temperature'
ParentCol = ParentFile+'_TotalColumn'
ParentWidth = ParentFile+'_WidthField'
ParentVel = ParentFile+'_VelocityField'

FileTemp = ParentTemp+'.FITS'
FileCol = ParentCol+'.FITS'
FileWidth = ParentWidth+'.FITS'
FileVel = ParentVel+'.FITS'

FileTempErr = ParentTemp+'_stderr.FITS'
FileColErr = ParentCol+'_stderr.FITS'
FileWidthErr = ParentWidth+'_stderr.FITS'
FileVelErr = ParentVel+'_stderr.FITS'

FITSparent = ''
ContFile = FITSparent+''

SaveData = True        # Save the data. This may overwrite the current data unless you change file names below.
SavePlots = True       # Save the plots as PNG images. 
SaveFITS = True        # Save the data as a FITS file

# Indicate the date and time of data adjustment.
ADJdate = date.today().strftime('%Y%m%d')
ADJtime = datetime.now().strftime('%H%M')

# Number of pixels
xpix = 513
ypix = 513

molectxt = '13CH3OH'


#####################
# SAVING PARAMETERS #
#####################

print 'Declaring saving parameters.'

if molectxt == '13CH3OH': molec = r'$^{13}$CH$_3$OH'
if molectxt == '13CH3CN': molec = r'$^{13}$CH$_3$CN'

FileTempOut = ParentTemp+'_ADJ_'+ADJdate+'_'+ADJtime+'.FITS'
FileColOut = ParentCol+'_ADJ_'+ADJdate+'_'+ADJtime+'.FITS'
FileWidthOut = ParentWidth+'_ADJ_'+ADJdate+'_'+ADJtime+'.FITS'
FileVelOut = ParentVel+'_ADJ_'+ADJdate+'_'+ADJtime+'.FITS'

FileTempErrOut = ParentTemp+'_stderr_ADJ_'+ADJdate+'_'+ADJtime+'.FITS'
FileColErrOut = ParentCol+'_stderr_ADJ_'+ADJdate+'_'+ADJtime+'.FITS'
FileWidthErrOut = ParentWidth+'_stderr_ADJ_'+ADJdate+'_'+ADJtime+'.FITS'
FileVelErrOut = ParentVel+'_stderr_ADJ_'+ADJdate+'_'+ADJtime+'.FITS'

FileTempPctErrOut = ParentTemp+'_pcterr_ADJ_'+ADJdate+'_'+ADJtime+'.FITS'
FileColPctErrOut = ParentCol+'_pcterr_ADJ_'+ADJdate+'_'+ADJtime+'.FITS'
FileWidthPctErrOut = ParentWidth+'_pcterr_ADJ_'+ADJdate+'_'+ADJtime+'.FITS'
FileVelPctErrOut = ParentVel+'_pcterr_ADJ_'+ADJdate+'_'+ADJtime+'.FITS'


def Load_FITS_File(FITSname):
    hdulist = fits.open(FITSname) # AstroPy FITS file loader, way faster than CASA, that steaming pile of junk, we should be loading the entire dataset into RAM
    data = hdulist[0].data[0] # Take only the Stokes I data; there isn't anything else anyway, but this simplifies things
    header = hdulist[0].header
    hdulist.close()
    data = np.swapaxes(data,0,2) # Switch axes so it is data[RA][DEC][spectrum]
    return data,header


"""
################
CODE STARTS HERE
################
"""

StartTime = time.time() # LETS GET READY TO RUUUMMMMMBLLE

print 'Loading original processed data for correcting.'
print ''

FITScol = Load_FITS_File(FileCol)
FITStemp = Load_FITS_File(FileTemp)
FITSwidth = Load_FITS_File(FileWidth)
FITSvel = Load_FITS_File(FileVel)

FITScolerr = Load_FITS_File(FileColErr)
FITStemperr = Load_FITS_File(FileTempErr)
FITSwidtherr = Load_FITS_File(FileWidthErr)
FITSvelerr = Load_FITS_File(FileVelErr)

FITScol = np.swapaxes(FITScol[0],0,2)
FITStemp = np.swapaxes(FITStemp[0],0,2)
FITSwidth = np.swapaxes(FITSwidth[0],0,2)
FITSvel = np.swapaxes(FITSvel[0],0,2)

FITScolerr = np.swapaxes(FITScolerr[0],0,2)
FITStemperr = np.swapaxes(FITStemperr[0],0,2)
FITSwidtherr = np.swapaxes(FITSwidtherr[0],0,2)
FITSvelerr = np.swapaxes(FITSvelerr[0],0,2)

FITScol = FITScol[0]
FITStemp = FITStemp[0]
FITSwidth = FITSwidth[0]
FITSvel = FITSvel[0]

FITScolerr = FITScolerr[0]
FITStemperr = FITStemperr[0]
FITSwidtherr = FITSwidtherr[0]
FITSvelerr = FITSvelerr[0]

ParName = ['TotalColumn','Temperature','WidthField','VelocityField']
ParFile = [FITScol,FITStemp,FITSwidth,FITSvel]
ParErrFile = [FITScolerr,FITStemperr,FITSwidtherr,FITSvelerr]
ParentFile = [ParentCol,ParentTemp,ParentWidth,ParentVel]
ParPrint = ['column','temperature','width','velocity']

#ParUnit = [r'cm$^{-2}$','K','MHz',r'km s$^{-1}$']
ParUnit = [r'cm$^{-2}$','K',r'km s$^{-1}$',r'km s$^{-1}$']


FileOut = [FileColOut,FileTempOut,FileWidthOut,FileVelOut]
FileErrOut = [FileColErrOut,FileTempErrOut,FileWidthErrOut,FileVelErrOut]
FilePctErrOut = [FileColPctErrOut,FileTempPctErrOut,FileWidthPctErrOut,FileVelPctErrOut]



discontpix = 0
for p in range(3):
    print 'Checking '+ParPrint[p]+' map for discontinuities caused by outlier pixels.'
    ppix=0
    for i in range(len(ParFile[p])):
        for j in range(len(ParFile[p][i])):
    
            areasum=0
            arean=0
            if (np.isfinite(ParFile[p][i][j])):
                for x in [(i-1),(i),(i+1)]:
                    for y in [(j-1),(j+1)]:
                        if np.isfinite(ParFile[p][x][y]):
                            areasum += ParFile[p][x][y]
                            arean += 1
                for x in [(i-1),(i+1)]:
                    for y in [j]:
                        if np.isfinite(ParFile[p][x][y]):
                            areasum += ParFile[p][x][y]
                            arean += 1
                if arean > 0:
                    areaave = areasum/arean
                    if p == 0 or p == 1: threshval = 0.10
		    if p == 2 or p == 3: threshval = 0.10
                    if (ParFile[p][i][j] < (1-threshval)*areaave) or (ParFile[p][i][j] > (1+threshval)*areaave):
                        discontpix += 1
                        ppix += 1
                        FITScol[i][j] = np.nan
                        FITScolerr[i][j] = np.nan
                        FITStemp[i][j] = np.nan
                        FITStemperr[i][j] = np.nan
                        FITSwidth[i][j] = np.nan
                        FITSwidtherr[i][j] = np.nan
                        FITSvel[i][j] = np.nan     
                        FITSvelerr[i][j] = np.nan
#    print str(ParPrint[p])+' = '+str(ppix)
print 'There were '+str(discontpix)+' outlier pixels removed for failing the discontinuity check.'
print ''





discontpix = 0
for p in range(len(ParFile)):
    print 'Checking '+ParPrint[p]+' map for pixels with high errors.'
    for i in range(len(ParFile[p])):
        for j in range(len(ParFile[p])):     
            nannum = 0
            if (np.isfinite(ParFile[p][i][j])) or (np.isfinite(ParErrFile[p][i][j])):
		if ((p == 0) and (np.abs(ParFile[p][i][j])<(3.0*np.abs(ParErrFile[p][i][j]))) or \
		   (p == 1) and (np.abs(ParFile[p][i][j])<(3.0*np.abs(ParErrFile[p][i][j]))) or \
		   (p == 2) and (np.abs(ParFile[p][i][j])<(3.0*np.abs(ParErrFile[p][i][j]))) or \
		   (p == 3) and (np.abs(ParFile[p][i][j])<(3.0*np.abs(ParErrFile[p][i][j])))):
			discontpix +=1
                        FITScol[i][j] = np.nan
                        FITScolerr[i][j] = np.nan
                        FITStemp[i][j] = np.nan
                        FITStemperr[i][j] = np.nan
                        FITSwidth[i][j] = np.nan
                        FITSwidtherr[i][j] = np.nan
                        FITSvel[i][j] = np.nan     
                        FITSvelerr[i][j] = np.nan
print 'There were '+str(discontpix)+' pixels removed for high errors.'         
print ''




"""
discontpix = 0
colpix = 0
temppix = 0
widthpix = 0
velpix = 0
for Index in range(len(ParFile)):
    print 'Manually removing visible outliers in '+ParPrint[p]+' map after previous discontinuity checks.'
    for i in range(len(ParFile[p])):
        for j in range(len(ParFile[p])):
            if ((Index == 1) and (ParFile[Index][i][j] < 50.0)): #\
               #((Index == 0) and (ParFile[p][i][j] > 9e16)) or \
               #((Index == 2) and (ParFile[Index][i][j] > 1.9)) or \
               #((Index == 3) and (ParFile[Index][i][j] > -2.9)):
                discontpix += 1
                if ((Index == 0)): colpix += 1
                if ((Index == 1)): temppix += 1
                if ((Index == 2)): widthpix += 1
                if ((Index == 3)): velpix += 1
                FITScol[i][j] = np.nan
                FITScolerr[i][j] = np.nan
                FITStemp[i][j] = np.nan
                FITStemperr[i][j] = np.nan
                FITSwidth[i][j] = np.nan
                FITSwidtherr[i][j] = np.nan
                FITSvel[i][j] = np.nan
                FITSvelerr[i][j] = np.nan
print 'There were '+str(discontpix)+' pixels manually removed as sources of discontinuity.'
print ' - Number of pixels manually removed from column density map: '+str(colpix)
print ' - Number of pixels manually removed from temperature map:    '+str(temppix)
print ' - Number of pixels manually removed from width map:          '+str(widthpix)
print ' - Number of pixels manually removed from velocity map:       '+str(velpix)
"""

discontpix = 0
for p in range(len(ParFile)):
    print 'Checking '+ParPrint[p]+' map for discontinuities caused by isolated pixels.'
    for i in range(len(ParFile[p])):
        for j in range(len(ParFile[p])):     
            nannum = 0
            if (np.isfinite(ParFile[p][i][j])):
                for x in [(i-1),(i),(i+1)]:
                    for y in [(j-1),(j),(j+1)]:
                        if np.isnan(ParFile[p][x][y]):
                            nannum += 1
                if nannum > 4:
                    discontpix += 1
                    FITScol[i][j] = np.nan
                    FITScolerr[i][j] = np.nan
                    FITStemp[i][j] = np.nan
                    FITStemperr[i][j] = np.nan
                    FITSwidth[i][j] = np.nan
                    FITSwidtherr[i][j] = np.nan
                    FITSvel[i][j] = np.nan
                    FITSvelerr[i][j] = np.nan

discontpix = 0
for i in range(len(ParFile[0])):
	for j in range(190):
	                discontpix += 1
                        FITScol[i][j] = np.nan
                        FITScolerr[i][j] = np.nan
                        FITStemp[i][j] = np.nan
                        FITStemperr[i][j] = np.nan
                        FITSwidth[i][j] = np.nan
                        FITSwidtherr[i][j] = np.nan
                        FITSvel[i][j] = np.nan     
                        FITSvelerr[i][j] = np.nan
print 'There were '+str(discontpix)+' outlier pixels removed for failing at noisy coordinates.'
print ''


print 'There were '+str(discontpix)+' isolated pixels removed for failing the discontinuity check.'         
print ''

print ''
print 'Converting velocity field from MHz to km/s units.'
print ''
FITSvel = -2.99*((FITSvel/1.559976025))
FITSvelerr = 2.99*((FITSvelerr/1.559976025))

print ''
print 'Converting width field from MHz to km/s units.'
print ''
FITSwidth = 2.99*((FITSwidth/1.559976025))
FITSwidtherr = 2.99*((FITSwidtherr/1.559976025))

ParFile = [FITScol,FITStemp,FITSwidth,FITSvel] # Redefining ParFile, just in case, to make sure adjustments carry over.
ParErrFile = [FITScolerr,FITStemperr,FITSwidtherr,FITSvelerr]
FITScol = np.ma.masked_where(np.isnan(FITScol),FITScol)
FITStemp = np.ma.masked_where(np.isnan(FITStemp),FITStemp)


print ''
print 'Plotting adjusted maps.'

for Index in range(len(ParName)):
    print ' - Plotting adjusted '+ParPrint[Index]+' map.'
    plt.figure()
    plt.title(molec+' '+ParPrint[Index])
    #plt.imshow(FITScol,cmap='rainbow')
    plt.pcolormesh(ParFile[Index],cmap='rainbow')
    cbar = plt.colorbar()
    cbar.set_label(ParUnit[Index])
    if (SavePlots): plt.savefig(ParentFile[Index]+'_ADJ.png')
    print ' - Plotting adjusted '+ParPrint[Index]+' error map.'
    plt.figure()
    plt.title(molec+' '+ParPrint[Index]+' error')
    plt.pcolormesh(ParErrFile[Index],cmap='rainbow')
    cbar = plt.colorbar()
    cbar.set_label(ParUnit[Index])
    if (SavePlots): plt.savefig(ParentFile[Index]+'_stderr_ADJ.png')
    print ' - Plotting adjusted '+ParPrint[Index]+' percent error map.'
    plt.figure()
    plt.title(molec+' '+ParPrint[Index]+' percent error')
    plt.pcolormesh(ParErrFile[Index]*100./ParFile[Index],cmap='rainbow')
    cbar = plt.colorbar()
    cbar.set_label('% error')
    if (SavePlots): plt.savefig(ParentFile[Index]+'_pcterr_ADJ.png')


if(SaveFITS):
    print 'Saving data as FITS files.'

    hdulist = fits.open(ContFile)
    header = hdulist[0].header
    hdulist.close()

    for Index in range(len(ParName)):
        print ' - Saving adjusted '+ParPrint[Index]+' data.'
        FakeMap = np.zeros((1,1,ypix,xpix))
        for i in range(ypix):
            for j in range(xpix):
                FakeMap[0][0][i][j] = ParFile[Index][i][j]
        hdu = fits.PrimaryHDU(FakeMap)
        hdu.header = header
        hdu.writeto(FileOut[Index])

        print' - Saving adjusted '+ParPrint[Index]+' error data.'
        FakeMap = np.zeros((1,1,ypix,xpix))
        for i in range(ypix):
            for j in range(xpix):
                FakeMap[0][0][i][j] = ParErrFile[Index][i][j]
        hdu = fits.PrimaryHDU(FakeMap)
        hdu.header = header
        hdu.writeto(FileErrOut[Index])

        print' - Saving adjusted '+ParPrint[Index]+' percent error data.'
        FakeMap = np.zeros((1,1,ypix,xpix))
        for i in range(ypix):
            for j in range(xpix):
                FakeMap[0][0][i][j] = ParErrFile[Index][i][j]*100./ParFile[Index][i][j]
        hdu = fits.PrimaryHDU(FakeMap)
        hdu.header = header
        hdu.writeto(FilePctErrOut[Index])

    

 
plt.show()
