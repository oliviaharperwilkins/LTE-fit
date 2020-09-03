# LTE-fit
The fitting.py script was written for mapping parameters such as rotational temperature, column density, velocity field, and line width from observations with the Atacama Large Millimeter/submillimeter Array (ALMA). The fitting consists of a least-squares fitting program called LMFIT and assumes local thermodynamic equilibrium (LTE). It is heavily based off of scripts by P. B. Carroll.

This file goes over the specifics of each section of the code, including the types of files needed and what the code is doing. Some short descriptions are included as comments in the actual code, but this README file offers a more comprehensive overview.

The sections of the code are as follows:
- Input files
- Run options
- Map parameters
- Save parameters
- Defining functions
- Fitting parameters
- Derive parameters

------------------------------------------------------------------------------------

INPUT FILES

This part of the code defines which files will be used for fitting. All image files should be in FITS formats. For CASA images, this can be achieved by using the CASA command 'exportfits'. 

The image cube files and catalogue files are listed in arrays FITSfilenames and MolDataFileNames. Even if only one file of a given type is used, it should be entered as a single entry in the array since the code is optimized to use arrays of file names. 
- Image cube files must be in FITS format.
- Catalogue files must be in .tsv format with the columns. These can be exported from Splatalogue.online and must have eight (8) columns for frequency (MHz), frequency uncertainty (MHz), CDMS/JPL intensity, Sij<mu>**2 (Debye**2), log(Einstein A coefficient), upper energy (K), "masK" (1 = include, 0 = exclude), and upper state degeneracy. (If using Splatalogue, be sure to remove any molecular formulas, molecule names, and quantum numbers.) Einstein A isn't needed, so if this number isn't available, you can just leave this blank.

Offsets (array): An array used for loading data sets of different sizes. Larger data sets (i.e. larger maps) will need to start shifted in. Each element in the array is the number of pixels for the corresponding FITS file (i.e. of the same index in FITSfilenames) to be shifted in. This can be circumvented (and all zeroes used in the array) but having data sets of the same size, for example by using the CASA command 'imsubimage'. 

MolDataTag (array): An array of integer tags to distinguish what molecules are being targeted in each spectral window. If a FITS image contains multiple lines to be modeled, list it twice in the FITSfilenames array, each of which have different corresponding MolDataTag entries.

MolDataTxt (array): An array of text strings showing the label for each MolDataTag. It is important that these strings are consistent for identical species. 

molec (string): Specifies the molecule of interest for titles of plots. To use TeX format, precede the string with the letter *r*.


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

------------------------------------------------------------------------------------

MAP PARAMETERS

Specify the pixels along the X and Y axes where to start and stop the mapping. This saves time by excluding regions where you already know there is no emission or emission interest. 

------------------------------------------------------------------------------------

SAVE PARAMETERS

FileDate and FileTime set the date and time at the initialization of the code.

ParentFile is the file name prefix for any output files saved during the running of the code. By default, file names are 'ALMA_fitting_<Date>_<Time>_<molectxt>_<parameter>'
  
  
------------------------------------------------------------------------------------

DEFINING FUNCTIONS

Partition_Function(T,molectxt): Returns the partition function pased on a power law fit of the partition function against temperature available in Splatalogue. The power is fixed to 1.5, which may not precisely match the fit of the actual data. Additional power law fits can be added as needed. 
- T = temperature
- molectxt = molecular species

LTE_Intensity(NT,Q,T,EUpper,v,Sij): Returns the integrated intensity for a line, based on equations in Remijan et al. (2005, ApJ, 626). This function is coded for single-dish emission and so is *not* used for ALMA data.
- NT = total column density
- Q = partition function
- T = temperature
- EUpper = upper energy
- v = frequency 
- Sij = transition intensity

LTE_Intensity_Interferometer: Returns the integrated intensity for a line, based on equations in Remijan et al. (2005, ApJ, 626). This function is coded for interferometers, such as ALMA, and assumes an input of Jy/beam (which is what ALMA uses for intensity). The beam argument should be supplied as (bmaj)(bmin); however, since the code extracts square regions, the beam size is simplified for a small CPU speed-up.
- NT = column density
- Q = partition function
- T = temperature
- EUpper = upper energy
- v = frequency
- Sij = transition intensity
- B = 0.5
- beam = beam size

MHz_to_Kms: Converts a frequency in MHz to velocity space in km/s.
- Frequency_MHz = frequency to be converted to velocity units
- CenterFrequency_MHz = center frequency of spectral window

Kms_to_MHz: Converts velocity in km/s to frequency in MHz.
- Kms = velocity to be converted to frequency units
- CenterFrequency_MHz = center frequency of spectral window

Load_Mol_Data: Loads in the catalogue data and puts it into a format useable for later in the code with six indeces: frequency, upper state energy, transition intensity, upper state degeneracy, Einstein A coefficient (not used in the LTE_Intensity_Interferometer function), and molecular species tag.
- Filename = tsv file for molecule
- Tag = integer tag for the individual molecule

Build_XAxis: Extracts the frequency axis information from the FITS header of the spectral window file.
- Header = header extracted from FITS file

LMFIT_residual: Residual procedure incorporated into optimizing LMFIT.
- params = physical parameters being fit
- x = frequency
- data = spectrum intensities
- CatalogueData = catalogue data
- NTCount = number of unique molecules for fits

Load_Data: Load data from fits file.
FileNames = array of FITS file names

Array_RMS: Get array of root-mean-squared uncertainty for second fitting procedure.

------------------------------------------------------------------------------------
FITTING PARAMETERS

This part of the code sets up the arrays that will hold physical parameter information, namely column density ("NT#"), temperature ("Temperature"), line width ("Sigma"), and velocity ("Xoffset"). At this point, the molecular line data and data image cubes are also loaded in. 

The data are combined into a single array. Arrays are set up to hold all X (frequency) arrays. A mask array is set up with all entries set to "False". For each line in the transition catalogue, frequencies close enough to the rest frequency, i.e. within the fitting width (FittingWidth_Freq), have the mask set to "True". Note that this can add issues if the limit is set near or lower than the doppler shift. Since the entire mask array starts off as having all False entries and is only adjusted to switch to True, lines close to one another will not be remasked if looped over again. Another array, FittingXAxis, is made to store the axes that will actually be used in the fit, based on the frequencies masked as true. The FinalXAxis is the FittingXAxis but as a numpy array. 

Using the AxisMask array, the FITS data file are looped over and combined into one array. Each band within the data is looped over and appended to the FinalFile, which is then turned into a numpy array. This is the sister array to FinalXAxis.

------------------------------------------------------------------------------------
DERIVE PARAMETERS

This part of the code starts by setting up a counter of pixels included in the fit (DetectedPixels) and sets the start time so the time to run the code can be calculated and printed at the conclusion of the script running. The initial processing starts by setting up data file strings and text files to which the derived parameters will be written.

The code goes through the image cube pixel-by-pixel and only attempts to fit pixels where the maximum intensity in the spectrum at that pixel exceeds the fit threshold set in [Run options] section of the script. The LMFIT function, a least-squares fit, is used to fit the data assuming LTE and is assumed to be a reasonable fit until shown otherwise (e.g. with error checks). Spectra for the designated test pixel and surrounding pixels are plotted and will be displayed at the conclusion of the code, showing the data with the LMFIT results plotted over the data. 

If ErrorCheck = True, pixels are thrown out if any of the following are true:
- The error is too high (i.e. the resulting value is less than three (3) times the standard error).
- The standard error is infinite.
- The standard error is 0.
- The standard error is not defined.
An optional full confidence interval calculation can be attempted as well, if Confidence = True.

Pixels that are successfully fitted and pass the error check are then transfered to ParameterMap and ErrorMap arrays that contain the four parameters at each pixel. The results are also written to a text file. After fits are attempted at all pixels, the axes are swapped to put the parameter index before the pixel indices. 

An optional discontinuity check is run and tosses out pixels if any of the following are true:
- The pixel is surrounded by too many NaNs.
- The pixel has a parameter value that is significantly higher or lower than the average of the surrounding pixels.
After this, all NaN pixels are masked and the completion time is recorded.

At the end of the script, the number of fittable pixels found are printed along with a summary of pixel rejections. Maps for parameters, if PlotPrimaryMaps = True, are generated and saved as PNG files if SavePlots = True. If PlotErrorMaps = True, maps of the standard error for each parameter are also plotted and saved as PNG files if SavePlots = True.

The resulting parameter maps can also be saved to FITS files if SaveFITS = True. FITS files are saved by extracting the header from the continuum FITS file to be applied to the parameter map FITS file. For each parameter, a temperary "FakeMap" is set up, and unmasked pixels in the parameter maps are transferred to the FakeMap, which is then written to a FITS file. The error data are also saved to FITS files in the same manner if SaveFITSstderr = True.

Plots, including the test spectra, primary maps (if selected), and error maps (if selected) are shown at the conclusion of the script. 
