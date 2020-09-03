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
SaveFITSerrBF = False   	# Save error data at best fit as FITS files. 
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



