"""
Created on Wed May 15 15:35:00 2019
@author: Olivia Harper Wilkins

Based on scripts by P. B. Carroll.

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

###############
# INPUT FILES #
###############


print ''
print 'Declaring input files.'

FitsParent = ''     # Path to directory containing FITS files containing spectral line data.
CatParent = ''      # Path to directory containing line catalogue files.

FITSfilenames = ['','',''] # Names of individual data files from spectral windows to be used.
Offsets = [0,0,0]
MolDataTag = [0,0,0] # e.g. 0 = 13CH3OH
MolDataTxt = ['13CH3OH','13CH3OH','13CH3OH']
molectxt = MolDataTxt[0]
molec = r'$^{13}$CH$_3$OH'
UniqueMolecules = len(set(MolDataTag)) # Convert the tag list to a set and count how long that set is.
print ' - There are '+str(UniqueMolecules)+' unique species being fitted.'

MolDataFileNames = ['','','']  # Names of catalogue files.
print 'Designating catalogue files.'

ContFile = '' # Continuum FITS image.

###############
# RUN OPTIONS #
###############

xpix = 513 # Number of pixels along X axis in subimage.
ypix = 513 # Number of pixels along Y axis in subimage.

print 'Loading run options.'

SaveData = True         	# Save the data. THIS WILL OVERWRITE THE CURRENT DATA UNLESS YOU CHANGE THE FILE NAMES BELOW.
ErrorCheck = True       	# Check and mask pixels with high error. 
SecondaryFitting = False	# Perform secondary fitting of the data.
SavePlots = True        	# Save the plots. This just saves literally anything you have enabled for plotting.
Confidence = False      	# Run LMFITs explicit confidence intervals.
SaveFITS = True         	# Save data as FITS files.
SaveFITSstderr = True   	# Save error data at 95% confidence as FITS files.
PlotPrimaryMaps = True  	# Plot NT and T_rot
PlotErrorMaps = True    	# Plot error maps
DiscontCheck = False    	# Discontinuity check for isolated pixels

itest = 277
jtest = 206

print 'Reading in fit settings.'

VLSR = 6.5                                   # Initial guess at VLSR
bmaj = 0.38                                  # [''] Size of the synthesized beam along the major axis.
bmin = 0.35                                  # [''] Size of the synthesized beam along the minor axis.
BeamSize = bmaj*bmin
CellSize = 0.05                              # [''] Size of a pixel.
BeamPixX = 7                                 # Pixels to use in X.
BeamPixY = 7                                 # Pixels to use in Y.
BeamSize = (BeamPixX*BeamPixY)*CellSize**2.0 # Superseding the true beam with the actual beam we use.
FitThreshold = 0.05                         # Threshold for initial pixel testing. 

WidthGuess = 0.04
XOffsetGuess = -3.3
TGuess = 100
NTGuess = 1.0E+15
WidthLimit = 2.0 # Only the upper limit is specified since, by definition, the width can never be below 0.
XOffsetLow = -5.3
XOffsetHigh = -1.3
TLow = 10.0
THigh = 400.0
NTLow = 5.0E+07
NTHigh = 1.0E+18

FittingWidth_Freq = 15.0
PointingCenter = SkyCoord('05h35m14.5s','-05d22m30.9s')

##################
# MAP PARAMETERS #
##################

StartX = 5             # Starting pixel in X.
StopX  = 510           # Stopping pixel in X.
StartY = 5 	           # Starting pixel in Y.
StopY  = 510           # Stopping pixel in Y.
SizeX  = StopX-StartX  # Total pixels in X.
SizeY  = StopY-StartY  # Total pixels in Y.


###################
# SAVE PARAMETERS #
###################

print 'Declaring saving parameters.'

FileDate = date.today().strftime('%Y%m%d')
FileTime = datetime.now().strftime('%H%M')

ParentFile = 'ALMA_fitting_'+FileDate+'_'+FileTime
FileName = ParentFile+'_'+molectxt+'_'


######################
# DEFINING FUNCTIONS #
######################

h = 6.626070040e-34 # [J s]  Planck constant
cms = 2.99792458e8  # [m/s]  speed of light
cc = 2.99792458e10  # [cm/s] speed of light
k = 1.38064852e-23  # [J/K]  Boltzmann constant

def Partition_Function(T,molectxt):
# Empirical partition function from power law fit of Splatalogue data. Power was fixed to 1.5.
# This may need to be adjusted for CH3OH and other species.
    if molectxt == '13CH3OH': # 13CH3OH
        return 0.39927*T**1.5
    if molectxt == 'CH3OD': 
        return 1.6277*T**1.5 
    if molectxt == 'CH2DOH: 
        return 0.5941*T**1.5 
    if molectxt == '13CH3CN':
        return 2.34051*T**1.5

def LTE_Intensity(NT,Q,T,EUpper,v,Sij):
# Compute an LTE intensity (technically the integrated intensity) for a line based on Remijan et al. (2005).
# Coded for single-dish emission.
    C = 4.8E-5
    A = (1-(exp((C*v)/T)-1)/(exp((C*v)/2.7)-1))
    return (NT*v*Sij*A)/((1.8E+14)*Q*exp(EUpper/T))

def LTE_Intensity_Interferometer(NT,Q,T,EUpper,v,Sij,B,beam):
# Compute an LTE intensity (technically the integrated intensity) for a line based on Remijan et al. (2005). The listed frequency has to be in GHz for this calculation, not MHz like it implies in the paper.
# Currently coded for an interferometer.
# This technically returns the integrated intensity of the line, assuming an input of Jy/beam (which is what ALMA uses).
# The beam argument should be supplied as bmaj*bmin. However, since the code currently extracts square regions, this should be exactly right.
# This simplification is done for a small cpu speed up.
    Term1 = (NT*B*(beam)/(2.04E+20))
    Term2 = ((((v/1000.0)**3.0)*Sij)/(Q*exp(EUpper/T)))
    return Term1*Term2

def MHz_to_Kms(Frequency_MHz,CenterFrequency_MHz):
# Converts frequency in MHz to velocity space in km/s
    return (Frequency_MHz/CenterFrequency_MHz)*2.99E+5

def Kms_to_MHz(Kms,CenterFrequency_MHz):
# Converts velocity space in km/s to frequency in MHz
    return (Kms/2.99E+5)*CenterFrequency_MHz

def Load_FITS_File(FITSname):
    hdulist = fits.open(FITSname) # AstroPy FITS file loader, way faster than CASA, that steaming pile of junk, we should be loading the entire dataset into RAM
    data = hdulist[0].data[0] # Take only the Stokes I data; there isn't anything else anyway, but this simplifies things
    header = hdulist[0].header
    hdulist.close()
    data = np.swapaxes(data,0,2) # Switch axes so it is data[RA][DEC][spectrum]
    return data,header


def Load_Mol_Data(Filename,Tag):
    Moldata = plb.loadtxt(Filename, skiprows = 1)
    Catalogue_Freq = Moldata[:,0]              # Frequencies
    Catalogue_Error = Moldata[:,1]             # Catalogue errors; we don't really use this.
    Catalogue_JPLStr = Moldata[:,2]            # JPL log(strength)
    Catalogue_SijMu = Moldata[:,3]             # Sij values
    Catalogue_A = Moldata[:,4]                 # Einstein A coefficient (not necessary)
    Catalogue_UpperState = Moldata[:,5]        # Upper state energies
    Catalogue_Mask = Moldata[:,6]              # To fit or not to fit
    Catalogue_Degeneracy = Moldata[:,7]        # Degeneracy; divide by this number for multiple lines to fit just a single line

    CatalogueLength = 0
    for i in range(len(Catalogue_Freq)):
        if (Catalogue_Mask[i] == 1):
            CatalogueLength += 1

    Catalogue = np.zeros((CatalogueLength,6)) # Catalogue[Frequency (0), Upper state energy (1), SijMu^2 (2) Degeneracy (3), Einstein A (5)][Line #]
    
    Count = 0
    for i in range(len(Catalogue_Freq)):
        if (Catalogue_Mask[i] == 1):
            Catalogue[Count][0] = Catalogue_Freq[i]
            Catalogue[Count][1] = Catalogue_UpperState[i]
            Catalogue[Count][2] = Catalogue_SijMu[i]
            Catalogue[Count][3] = Catalogue_Degeneracy[i]
            Catalogue[Count][4] = 10.0**Catalogue_A[i]
            Catalogue[Count][5] = Tag
            Count += 1

    return Catalogue

def Build_XAxis(Header):
    if (Header['CTYPE3'] != 'FREQ'):
        return -1
    StartX = Header['CRVAL3']
    XStep = Header['CDELT3']
    X = np.zeros(int(Header['NAXIS3']))
    for i in range(int(Header['NAXIS3'])):
        X[i] = (StartX+i*XStep)/1.0E+6
    return X


def LMFIT_residual(params,x,data,CatalogueData,NTCount): 
    y = np.zeros_like(x)
    NT = []
    for i in range(NTCount):
        NT.append(params['NT%d' % i])
    sigma = params['Sigma']
    xoffset = params['Xoffset']
    T = params['Temperature']
    BaseLine = params['BaseLine']
    LineCount = 0
    for i in range(NTCount):
        Line = CatalogueData[LineCount]
        while ((Line[5] == i) and (LineCount < len(CatalogueData))):
            Line = CatalogueData[LineCount]
            Partition = Partition_Function(T,molectxt)
#            y = y + (1.0/(sigma*2.35482))*LTE_Intensity_Interferometer(NT[i],Partition,T,Line[1],Line[0],Line[2],0.5,BeamSize)*exp(-(x-Line[0]-xoffset)**2.0/(2.0*sigma**2.0))
            y = y + LTE_Intensity_Interferometer(NT[i],Partition,T,Line[1],Line[0],Line[2],0.5,BeamSize)*exp(-(x-Line[0]-xoffset)**2.0/(2.0*sigma**2.0))
            LineCount += 1
    y += BaseLine
    return (data-y)


def Load_Data(FileNames,StartX,StopX,StartY,StopY,Offsets):
    FITSfiles = []
    FITSfiles.append([])
    FITSfiles.append([])
    for Index,Name in enumerate(FileNames):
        Data,Header = Load_FITS_File(Name)
        FITSfiles[0].append(np.copy(Data[StartX+Offsets[Index]:StopX+Offsets[Index],StartY+Offsets[Index]:StopY+Offsets[Index]]))
        FITSfiles[1].append(Header)
        del Data,Header
    FITSfiles = np.array(FITSfiles)
    return FITSfiles[0],FITSfiles[1]

def Array_RMS(Array):
    Temp = Array**2.0
    Avg = np.average(Temp)
    return np.sqrt(Avg)
  
"""
CODE STARTS HERE  
"""
StartTime = time.time()

######################
# FITTING PARAMETERS #
######################

print ''
print 'Setting fitting parameters.'

params = Parameters()
for i in range(UniqueMolecules):
    params.add("NT%d" % i, value = NTGuess, min = NTLow, max = NTHigh)
params.add("Sigma", value = WidthGuess, min = 0.0, max = WidthLimit)
params.add("Xoffset", value = XOffsetGuess, min = XOffsetLow, max = XOffsetHigh)
params.add("Temperature", value = TGuess, min = TLow, max = THigh)
params.add("BaseLine", value = 0.0, min = -10.0, max = 10.0, vary = False)
ParameterKey = []

for i in range(UniqueMolecules):
    ParameterKey.append("NT%d" % i)
ParameterKey.append("Sigma")
ParameterKey.append("Xoffset")
ParameterKey.append("Temperature")
ParameterKey.append("BaseLine")

print ''
print 'Loading molecular line data.'

MolData = []
for Index,FileName in enumerate(MolDataFileNames): 
    Temp = Load_Mol_Data(FileName,MolDataTag[Index])
    for Line in Temp: MolData.append(Line)
MolData = np.array(MolData)
print 'Loaded %d molecular lines.' % len(MolData)

print ''
print 'Loading FITS files.'
FITSdata,FITSheaders = Load_Data(FITSfilenames,StartX,StopX,StartY,StopY,Offsets) # Data are trimmed in RA/DEC (really in pixel space) to minimize RAM usage.
print 'Loaded %d FITS files.' % len(FITSdata)
for Band in FITSdata:
    for i in range(SizeX):
        for j in range(SizeY):
            Band[i][j] -= np.average(Band[i][j])*0.5#0.75
if molectxt == '13CH3OH':
    for i in range(SizeX):
        for j in range(SizeY):
            FITSdata[1][i][j] += 0#0.00186
    for i in range(SizeX):
        for j in range(SizeY):
            #        FITSdata[0][i][j] += 0.006
            FITSdata[1][i][j] += 0.004
            FITSdata[2][i][j] += 0.0005
if molectxt == '13CH3CN':
    for i in range(SizeX):
        for j in range(SizeY):
            FITSdata[0][i][j] += 0.004

print ''
print 'Combining data into a single object.'

XAxes = [] # An array to hold all of our X arrays.
for Header in FITSheaders:
    XAxes.append(Build_XAxis(Header))

AxisMask = [] # A mask array to track the parts of the array.
for Array in XAxes:
    AxisMask.append(np.zeros_like(Array,dtype=bool)) # Built a second list of numpy boolean arrays initialized to False, since it is an array of 0s which is False when it is specified as a boolean.

nTrue=0
for Line in MolData: # Loop over lines in the full molecular line data set.
    LineFrequency = Line[0] # Set the line frequency as the frequency stored in the line object. 
    for BandIndex, Band in enumerate(XAxes): # Loop over each individual X axis
        for PointIndex,Point in enumerate(Band): # Run through each point within an individual axis. This is just a more user-friendly, though inefficient, way of looping over all the relevant slices.
            if (np.abs(LineFrequency - Point) < FittingWidth_Freq) and (AxisMask[BandIndex][PointIndex] == False) and (PointIndex < len(Band)-1): 
                AxisMask[BandIndex][PointIndex] = True
                nTrue+=1

FittingXAxis = [] # New array for storing the axes we'll actually use to do the fits. 
for BandIndex,Band in enumerate(AxisMask): # Loop over the mask, which is the same as looping over the original axes.
    for PointIndex,Point in enumerate(Band): # Loop over individual points/slices within the band.
        if (Point == True):
            FittingXAxis.append(XAxes[BandIndex][PointIndex])

FinalXAxis = np.asarray(FittingXAxis) 
del FittingXAxis 
        
DataSlicing = [] # This will eventually look like DataSlicing[Band][Window][A tuple that has the Start and Stop point of each window]
for BandIndex,Band in enumerate(AxisMask): # Loop over each FITS file's X axis
    Count = 0
    DataSlicing.append([]) # Make a list of start/stop points for each axis.
    while (Count < len(Band) - 1): 
        while((Count < len(Band) - 1) and (Band[Count] == False )):
            Count += 1 #ONWARD!
        if (Count == len(Band) - 1): continue # If this takes us to the end of the array, break out of the loop and let the next 'while' loop take us into the next band without adding a start/stop point to our list, basically we found nothing more to use
        Start = Count
        while ((Band[Count] == True) and (Count < len(Band) - 1)): # Now loop through the group of TRUEs.
               Count += 1
             #  print Count
        if Count == len(Band)-1:
            Stop = Count
        else:
            Stop = Count
#        print 'Start = '+str(Start)
#        print 'Stop = '+str(Stop)
        DataSlicing[BandIndex].append([Start,Stop]) # Add this start/stop tuple to our list and then head back to the top to loop through the next region of Falses.

FinalFile = [] 
XPoints = 0
for BandIndex,Band in enumerate(DataSlicing):
    for WindowIndex,Window in enumerate(Band):
        XPoints += np.abs(Window[1]-Window[0])
        for Slice in range(Window[0],Window[1]):
            FinalFile.append(FITSdata[BandIndex][:,:,Slice]) .
FinalFile = np.asarray(FinalFile) 
FinalFile = np.swapaxes(FinalFile,0,2)
FinalFile = np.swapaxes(FinalFile,0,1)
del FITSdata

print ''
print '------------------------------'

#####################
# DERIVE PARAMETERS #
#####################

DetectedPixels = 0 # A count of pixels where we decided to fit.
StartTime = time.time()
print ''
print 'Starting initial processing.'
print ''
FileString = ''
FileHandle = open("%s.log" % ParentFile, 'w')
FileString += "X:\tY:\t"
for Par in ParameterKey: FileString += "%s\t%s Err\t" % (Par,Par)
FileString += "%s\tNFEV\t" % ("Chi Squared")
FileString += "\n"
FileHandle.write(FileString)
FileHandle.close()

DataFileString = ''
DataFileHandle = open("%s_StdErr.log" % ParentFile, 'w')
DataFileString+= "X:\tY:\t"
for Par in ParameterKey: DataFileString += "%s\t%s Err\t" % (Par,Par)
DataFileString += "\n"
DataFileHandle.write(DataFileString)
DataFileHandle.close()

errhigh = 0
errnan  = 0
errzero = 0
errundef = 0
errci = 0

ParName = ['TotalColumn','Temperature','WidthField','VelocityField']
ParAlias = ['NT0','Temperature','Sigma','Xoffset'] # Alias used in LMfit parameterization script.
ParPrint = ['column','temperature','width','velocity']
ParUnits = [r'cm$^{-2}$','K','MHz','km/s']
ParPrintPlot = ['Total Column','Temperature','Width Field','Velocity Field']
ParameterMap = np.zeros((ypix,xpix,len(ParName)))
ParameterMap[:][:][:] = np.nan
ErrorMap = np.zeros((ypix,xpix,len(ParName)))
ErrorMap[:][:][:] = np.nan




for i in range(SizeX):
    for j in range(SizeY):
        print 'Processing pixel X: %d Y: %d' % (i+StartX,j+StartY)
        if (np.amax(FinalFile[i][j]) > FitThreshold):
            DetectedPixels += 1 # Track useful pixels.
            Mini = lmfit.Minimizer(LMFIT_residual, params, fcn_args = (FinalXAxis,FinalFile[i][j],MolData,UniqueMolecules))
            Result = Mini.minimize(maxfev = 500)
            if (SecondaryFitting): # If we want to double-check our work, we can set this to TRUE to do additional fit checks.
                Result2 = Mini.minimize(method = "basinhopping") # Currently using the basinhopping algorithm. WARNING: this thing takes 100-1000x more function evaulations per fit and is really time consuming. It also does a really nice job of exploring phase space.
                if (Array_RMS(Result2.residual) < Array_RMS(Result.residual)): # Take the result with the best residual since basinhopping doesn't seem to derive a proper Chi-Squared.
                    Result = Result2
                    for Index,Par in enumerate(params): params[Par].value = Result.params[Par].value
                    Result = Mini.minimize(maxfev = 500)
                    for i in range(UniqueMolecules):
                        params["NT%d" % i].value = NTGuess
                    params["Sigma"].value = WidthGuess
                    params["Xoffset"].value = XOffsetGuess
                    params["Temperature"].value = TGuess
                    params["BaseLine"].value = 0.0
            GoodFit = True # Innocent until proven guilty.

            for tx in [itest, itest-1, itest-3, itest+4]:
                for ty in [jtest, jtest+1, jtest+3, jtest-4]:
                    if i == (tx-StartX) and (j == ty-StartY):
                        xplot = np.arange(min(FinalXAxis),max(FinalXAxis),0.2)
                        yplot = np.zeros_like(xplot)
                        for m in range(0,len(MolData)):
                            yplot = yplot + LTE_Intensity_Interferometer(Result.params['NT0'].value,Partition_Function(Result.params['Temperature'].value,molectxt),Result.params['Temperature'].value,MolData[m][1],MolData[m][0],MolData[m][2],0.5,BeamSize)*exp(-(xplot-MolData[m][0]-Result.params['Xoffset'].value)**2/(2*Result.params['Sigma'].value**2))
                        print 'Plotting test spectrum for pixel X: '+str(itest)+' Y: '+str(jtest)
                        print 'NT = '+str("{:.2e}".format(Result.params['NT0'].value))+' cm**-2'
                        print 'Trot = '+str(Result.params['Temperature'].value)+' K'
                        plt.figure()
                        plt.title('Test spectrum for pixel X: '+str(tx)+' Y: '+str(ty))
                        plt.ylabel('Relative intensity')
                        plt.xlabel('Frequency [MHz]')
                        ax = plt.gca()
                        ax.get_xaxis().get_major_formatter().set_useOffset(False)
                        plt.plot(FinalXAxis,FinalFile[i][j])
                        plt.plot(xplot,yplot,'r--')


            if (ErrorCheck):
                try:
                    for Index,Par in enumerate(Result.params):
                        if (Result.params[Par].vary):
                            if (np.abs(Result.params[Par].value) < (3.0*Result.params[Par].stderr)):
                                GoodFit = False
                                print ' %s error too high.' % Par
                                errhigh += 1
                            if (not np.isfinite(Result.params[Par].stderr)):
                                GoodFit = False
                                print ' %s error is non-finite.' % Par
                                errnan += 1
                            if (Result.params[Par].stderr < 1.0e-10):
                                GoodFit = False
                                print ' %s error is zero.' % Par
                                errzero += 1

                except: 
                     GoodFit = False
                     print ' %s error is not defined.' % Par
                     errundef += 1

            if (Confidence and Result.success and GoodFit): # Perform a full confidence interval calculation.
                try: ci = lmfit.conf_interval(Mini,Result) # Confidence interval
                except: 
                    GoodFit = False
                    print 'Could not calculate full confidence interval.'
                    errci += 1
                
            print Result.nfev,GoodFit
            if (GoodFit == True):
                


                ParameterMap[i+StartX][j+StartY][0] = Result.params['NT0'].value
                ParameterMap[i+StartX][j+StartY][1] = Result.params['Temperature'].value
                ParameterMap[i+StartX][j+StartY][2] = Result.params['Sigma'].value
                ParameterMap[i+StartX][j+StartY][3] = Result.params['Xoffset'].value

                ErrorMap[i+StartX][j+StartY][0] = Result.params['NT0'].stderr
                ErrorMap[i+StartX][j+StartY][1] = Result.params['Temperature'].stderr
                ErrorMap[i+StartX][j+StartY][2] = Result.params['Sigma'].stderr
                ErrorMap[i+StartX][j+StartY][3] = Result.params['Xoffset'].stderr
 
		FileString = ''
		FileHandle = open('%s_StdErr.log' % ParentFile,'a')
		FileString += '%d\t$d\t' % (i+StartX,j+StartY)
		for Par in ParameterKey:
			if (Par[0] == "N"):
				FileString += "%.3e\t%.3e\t" % (Result.params[Par].value, Result.params[Par].stderr)
			else:
				if (Par[0] != "B"):
					FileString += "%.3f\t%.3f\t" % (Result.params[Par].value,Result.params[Par].stderr)
		FileHandle.write(FileString)
		FileHandle.close()

                if (Confidence):
                    FileString = ''
                    FileHandle = open('%s.log' % ParentFile,'a')
                    FileString += '%d\t%d\t' % (i+StartX,j+StartY)
                    for Par in ParameterKey:
                        if (Par[0] == "N"):
                            FileString += "%.3e\t%.3e\t%.3e\t" % (Result.params[Par].value,ci[Par][2][1],ci[Par][4][1])
                        else:
                            if (Par[0] != "B"):
                                FileString += "%.3f\t%.3f\t%.3f\t" % (Result.params[Par].value,ci[Par][2][1],ci[Par][4][1])
                    FileString += "%.4f\t" % (Result.chisqr)
                    FileString += "%d" % (Result.nfev)
                    FileString += "\n"
                    FileHandle.write(FileString)
                    FileHandle.close()
                else: 
                    FileString = ''
                    FileHandle = open("%s.log" % ParentFile,'a')
                    FileString += "%d\t%d\t" % (i+StartX,j+StartY)
                    for Par in ParameterKey:
                        if (Par[0] == "N"):
                            FileString += "%.3e\t%.3e \t" % (Result.params[Par].value,Result.params[Par].stderr)
                        else:
                            FileString += "%.3e\t%.3e \t" % (Result.params[Par].value,Result.params[Par].stderr)
                    FileString += "%.4f\t" % (Result.chisqr)
                    FielString = "%d" % (Result.nfev)
                    FileString += "\n"
                    FileHandle.write(FileString)
                    FileHandle.close()                   

            else:
                for Index in range(len(ParAlias)):
                    ParameterMap[i+StartX][j+StartY][Index] = np.nan
                    ErrorMap[i+StartX][j+StartY][Index] = np.nan


ParameterMap = np.swapaxes(ParameterMap,0,2)
ErrorMap = np.swapaxes(ErrorMap,0,2)

if (DiscontCheck):
    print ''
    print 'Checking parameter maps for discontinuities in pixel profiles.'
    print ''
    discontpix=0
    for i in range(SizeX):
        for j in range(SizeY):     
            nannum = 0
            if (np.isfinite(ParameterMap[0][i+StartX][j+StartY])):
                for x in [(i-1),(i),(i+1)]:
                    for y in [(j-1),(j),(j+1)]:
                        if np.isnan(ParameterMap[0][x+StartX][y+StartY]):
                            nannum += 1
                if nannum > 4:
                    discontpix += 1
                    for p in range(4):
                        ParameterMap[p][i+StartX][j+StartY] = np.nan

    for i in range(SizeX):
        for j in range(SizeY):
            areasum = 0
            arean = 0
            if (np.isfinite(ParameterMap[1][i+StartX][j+StartY])):
                for x in [(i-1),(i),(i+1)]:
                    for y in [(j-1),(j+1)]:
                        if np.isfinite(ParameterMap[1][x+StartX][y+StartY]):
                            areasum += ParameterMap[1][x+StartX][y+StartY]
                            arean += 1
                for x in [(i-1),(i+1)]:
                    if np.isfinite(ParameterMap[1][x+StartX][j+StartY]):
                        areasum += ParameterMap[1][x+StartX][y+StartY]
                        arean += 1
                if arean > 0:
                    areaave = areasum/arean
                    threshval = 0.5
                    if (ParameterMap[1][i+StartX][j+StartY] > (1+threshval*areaave)) or ((1/(1-threshval))*ParameterMap[1][i+StartX][j+StartY] < areaave):
                        discontpix += 1
                        for p in range(4):
                            ParameterMap[p][i+StartX][j+StartY] = np.nan
                if arean == 0: 
                    for p in range(4):
                        ParameterMap[p][i+StartX][j+StartY] = np.nan
                                        
    print "A total of %d values failed the discontinuity check." % (discontpix) 

ParameterMap = np.ma.masked_where(np.isnan(ParameterMap),ParameterMap)

FullTime = time.time()
print ''
print '------------------------------'
print ''
print "There were %d (%.1f%%) fittable pixels found." % (DetectedPixels, (100.0*DetectedPixels)/(SizeX*SizeY))
if (ErrorCheck):
    print "A total of %d values failed the error check." % (errhigh+errnan+errzero)
    print "     - The error was too high for %d values." % errhigh
    print "     - The error was NaN for %d values." % errnan
    print "     - The error was 0 for %d values." %errzero
    print "     - The error was undefined for %d values." %errundef
    print "     - The full confidence interval could not be calculated for %d values." %errci
    print ""
#    print 'This leaves %d pixels that passed the error check.' % (DetectedPixels-(errhigh+errnan+errzero))
print ''
print "All processing complete."
print "Full execution time: "+str(FullTime-StartTime)+" s ("+str((FullTime-StartTime)/60.0)+" min)"

FileString = ""
FileHandle = open("%s.log" % ParentFile, 'a')
FileString += "\n"
FileString += "There were %d (%.1f%%) fittable pixels found." % (DetectedPixels, (100.0*DetectedPixels)/(SizeX*SizeY))
FileString += ''
FileString += "All processing complete."
FileString += "Full execution time: "+str(FullTime-StartTime)+" s ("+str((FullTime-StartTime)/60.0)+" min)"
FileHandle.write(FileString)
FileHandle.close()

if(PlotPrimaryMaps):
    print ''
    print 'Plotting primary maps.'
    
    for Index in range(len(ParName)):
        print ' - Plotting '+ParPrint[Index]+' map.'
        plt.figure()
        plt.title(molec+' '+ParPrintPlot[Index])
        plt.pcolormesh(ParameterMap[Index],cmap='rainbow')
        cbar = plt.colorbar()
        cbar.set_label(ParUnits[Index])
        if (SavePlots): plt.savefig(FileName+'_'+ParName[Index]+'Map.png')


if(PlotErrorMaps):
    print ''
    print 'Plotting error maps.'
    
    for Index in range(len(ParName)):
        print ' - Plotting '+ParPrint[Index]+' error map.'
        plt.figure()
        plt.title(molec+' '+ParPrintPlot[Index]+' Error')
        plt.pcolormesh(ErrorMap[Index],cmap='rainbow')
        cbar = plt.colorbar()
        cbar.set_label(ParUnits[Index])
        if (SavePlots): plt.savefig(FileName+'_'+ParName[Index]+'Map_Error.png')

if(SaveFITS):
    print ''
    print 'Saving data as FITS files.'
    #
    hdulist = fits.open(ContFile)
    header = hdulist[0].header
    hdulist.close()
    #
    for Index in range(len(ParName)):
        print ' - Saving '+ParPrint[Index]+' data.'
        FakeMap = np.zeros((1,1,ypix,xpix))
        for i in range(ypix):
            for j in range(xpix):
                if ((i >= StartY) and (i < StopY) and (j >= StartX) and (j < StopX)):
                    FakeMap[0][0][i][j] = ParameterMap[Index][i][j]
                # else: 
                #     FakeMap[0][0][i][j] = np.nan
        hdu = fits.PrimaryHDU(FakeMap)
        hdu.header = header
        hdu.writeto('plotdata/'+ParentFile+'_'+ParName[Index]+'.FITS')
    
if (SaveFITSstderr):
    print ''
    print 'Saving standard errors as FITS files.'

    ErrorName = ['TotalColumn_stderr','Temperature_stderr','WidthField_stderr','VelocityField_stderr']

    hdulist = fits.open(ContFile)
    header = hdulist[0].header
    hdulist.close()

    for Index in range(len(ErrorName)):
        print ' - Saving '+ParPrint[Index]+' error data.'
        FakeMap = np.zeros((1,1,ypix,xpix))
        for i in range(ypix):
            for j in range(xpix):
                if ((i >= StartY) and (i < StopY) and (j >= StartX) and (j < StopX)):
                    FakeMap[0][0][i][j] = ErrorMap[Index][i][j]
                else: FakeMap[0][0][i][j] = np.nan
        hdu = fits.PrimaryHDU(FakeMap)
        hdu.header = header
        hdu.writeto('plotdata/'+ParentFile+'_'+ErrorName[Index]+'.FITS')


plt.show()


