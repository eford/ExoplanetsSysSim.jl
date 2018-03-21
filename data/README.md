# This directory contains data files likely to be commonly used with SysSim
For more info about the file contents, see:  
http://exoplanetarchive.ipac.caltech.edu/docs/Kepler_stellar_docs.html  
http://exoplanetarchive.ipac.caltech.edu/docs/API_keplerstellar_columns.html

# Downloading stellar catalog files
Some csv files were too big for github, so you should download them yourself.  
The commands in download_stellar_tables.sh were used to download the table of stellar properties in CSV format.  
Reading CSV files is significantly slower than JLD files, so it's recommended to make a jld version for files to be used often.

# Downloading koi catalog files
We include the csv version that we've downloaded in the ExoplanetSySim repository  
I don't know how to script the downloading of the KOI tables from:  
http://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=koi  
When you download via the interactive table, I suggest selecting all rows and all columns to make sure you don't miss something  
See info at:  
http://exoplanetarchive.ipac.caltech.edu/docs/API_kepcandidate_columns.html  
http://exoplanetarchive.ipac.caltech.edu/docs/Kepler_KOI_docs.html  
http://exoplanetarchive.ipac.caltech.edu/docs/PurposeOfKOITable.html

Q: What are the contents of the JLD files?  
A: The JLD files are repackaged stellar catalog data from the different csv versions available on the Exoplanet Archive.  
   Some JLD files have been filtered and labeled with the appropriate suffix:  
* christiansen = Filtered on the FGK selection in the Christiansen et al. (2015) study [Teff = 4000-7000K; log g > 4; observed at least once in Q1-Q12 (only applies to Q1-Q16 file)]
* cks = Incorporates improved stellar parameters from the California Kepler Survey
* stellar_type = Filtered on stellar type as defined by the SAG13 working group [K: 3900K <= Teff < 5300K; G: 5300K <= Teff < 6000K; F: 6000 <= Teff < 7300K]

Q: What is the difference between q1_q16_koi_cand.csv and q1_q16_koi.csv?  
A: Danley Hsu generated q1_q16_koi_cand.csv based on q1_q16_koi.csv as well as previous KOI catalogs.  The primary difference is that he has used disposition information from previous KOI catalogs that was not reflected in the raw Q16 KOI catalog to select all KOIs with the "candidate" disposition.
