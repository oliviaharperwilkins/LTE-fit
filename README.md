# LTE-fit
The fitting.py script was written for mapping parameters such as rotational temperature, column density, velocity field, and line width from observations with the Atacama Large Millimeter/submillimeter Array (ALMA). The fitting consists of a least-squares fitting program called LMFIT and assumes local thermodynamic equilibrium (LTE).

This file goes over the specifics of each section of the code, including the types of files needed and what the code is doing. Some short descriptions are included as comments in the actual code, but this README file offers a more comprehensive overview.

The sections of the code are as follows:
- Input files
- Run options



------------------------------------------------------------------------------------

INPUT FILES

This part of the code defines which files will be used for fitting. All image files should be in FITS formats. For CASA images, this can be achieved by using the CASA command 'exportfits'. 

The image cube files and catalogue files are listed in arrays FITSfilenames and MolDataFileNames. Even if only one file of a given type is used, it should be entered as a single entry in the array since the code is optimized to use arrays of file names. 
- Image cube files must be in FITS format.
- Catalogue files must be in .tsv format with the columns. These can be exported from Splatalogue.online and must have eight (8) columns for frequency (MHz), frequency uncertainty (MHz), CDMS/JPL intensity, Sij<mu>**2 (Debye**2), log(Einstein A coefficient), upper energy (K), "masK" (1 = include, 0 = exclude), and upper state degeneracy. (If using Splatalogue, be sure to remove any molecular formulas, molecule names, and quantum numbers.) Einstein A isn't needed, so if this number isn't available, you can just leave this blank.

Offsets (array): An array used for loading data sets of different sizes. Larger data sets (i.e. larger maps) will need to start shifted in. Each element in the array is the number of pixels for the corresponding FITS file (i.e. of the same index in FITSfilenames) to be shifted in. This can be circumvented (and all zeroes used in the array) but having data sets of the same size, for example by using the CASA command 'imsubimage'. 

MolDataTag (array): An array of integer tags to distinguish what molecules are being targeted in each spectral window. If a FITS image contains multiple lines to be modeled, list it twice in the FITSfilenames array, each of which have different corresponding MolDataTag entries.

MolDataTxt (array): An array of text strings showing the label for each MolDataTag. It is important that these strings are consistent for identical species. 


------------------------------------------------------------------------------------

RUN OPTIONS

xpix, ypix (integers): The number of pixels along the X and Y axes in the (sub)image used in the line fitting.

There is a list of run options that you can designate as either True or False:
SaveData: if True, this saves the data and overwrites the current data set UNLESS the file names deeper in the code are changed.
ErrorCheck: if True, the script will check and mask pixels with high errors
SecondaryFitting: if True, this performs a secondary fitting of the data *NB that I have not used this so don't know whether it works* 
Confidence: if True, runs LMFIT's explicit confidence intervals.
SaveFITS: if True, saves the fitted parameter maps to FITS files.
SaveFITSstderr: if True, saves the error data at 95% confidence to FITS files.
SaveFITSerrBF: if True, saves teh error data at best fit to FITS files.
PlotPrimaryMaps: if True, plots column density and rotational temperature maps.
PlotErrorMaps: if True, plots error maps.
DiscontCheck: if True, runs a discontinuity check for isolated pixels. *NB that this is not robust and so there is a separate script for adjusting the results to account for isolated pixels or pixels with high uncertainties.*

itest,jtest (integers): Designates pixels for plotting test spectra to see whether the fits seem reasonable.

VLSR (float): Initial guess at the local standard of rest (LSR) velocity in km/s.

bmaj,bmin(floats): Size of the synthesized beam along the major and minor axes in arcseconds.

CellSize (float): Size of a pixel in arcseconds.

BeamPixX, BeamPixY (integers): Number of pixels to use along the X and Y axes. This should always be an odd number to place a single pixel at the center. This allows us to go pixel by pixel more smoothly by considering the surrounding pixels in the beam.

FitThreshold (float): This is the threshold for what pixels will have spectral line fits attempted. If the entire spectrum at a given pixel is below this threshold, it will be ignored. This saves time by not attempting to fit pixels where there is negligible emission.

WidthGuess (float): Initial guess of FWHM line width in MHz. Must be positive.

XOffsetGuess (float): Initial guess of the frequency shift (in MHz) of spectral lines. Can be positive or negative.

TGuess (float): Initial guess of the temperature in Kelvin.

NTGuess (float): Initial guess of the total column density for a molecule.

WidthLimit (float): Upper limit of line width to minimize the chances of the script fitting noise or neighboring lines. Since the width can never be below 0, there is no need to specify a lower limit width.

XOffsetLow,XOffsetHigh (floats): Allowed ranges of frequency shift of spectral lines, again to minimize the chances of the script fitting noise or neighboring lines.

TLow, THigh (floats): Allowed ranges of temperature for relevant interstellar environments to prevent the script from suggesting outlandish temperatures resulting from poor fits.

NTLow,NTHigh (floats): Allowed ranges of total column density for relevant interstellar environments to prevent the script from suggesting outlandish abundances resulting from poor fits.

FittingWidth_Freq (float): Bandwidth in MHz for test spectrum plots to make it easier to see the respective lines and surrounding noise and neighboring lines.

PointingCenter: The right ascension and declination used as the pointing center in the observations.
