The adjust_fit.py script is used to run error and discontinuity checks on the products from the fitting.py script. The error and discontinuity checks in the fitting.py script are not all that robust, so it is advisable to run the adjust_fit.py script, perhaps several times to try different thresholds.

This script consists of the following sections:
- Inputs
- Saving parameters

-------------------------------------------------------------------
INPUTS

Parent Dir (string): Path to directory containing parameter map files.

ParentDate (integer): Date on which parameter map generation was initialized in YYYYMMDD format, e.g. ParentDate = '20200516'.

ParentTime (integer): Time at which parameter map generation was initialized in HHMM format e.g. ParentTime = '2138'

FITSparent (string): Path to directory containing continuum file.

The FileX inputs probably do not need to be changed unless you have different prefixes or suffixes in your data files. If not, don't check these.

ContFile (string): Continuum file, which must be a FITS file.

There are three saving options that can be set to True:
- SaveData: if True, the adjusted data is saved to text files.
- SavePlots: if True, plots of the adjusted data are saved as PNG images.
- SaveFITS: if True, the adjusted data is saved to FITS files.

xpix,ypix (integers): Number of pixels in the X and Y dimensions.

molectxt (string): Molecular formula in plain format.

-------------------------------------------------------------------
SAVING PARAMETERS

Set up if-then statements for designating the TeX rich-text chemical formula. If molectxt == '<plain formula>': molec = r'<TeX formula>'.
