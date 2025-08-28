This is a code repository for the manuscript “Estimates of under-five mortality from a mobile phone survey during the COVID-19 outbreak: an observational study comparing three instruments in Malawi”, which has been accepted for publication in Demography.

The MATLAB file RaMMPS_U5M.m produces all figures and tables included in the manuscript, using some input data files: RaMMPS.csv, RaMMPSdst.csv, RaMMPSHH.csv, MICSmalawi.csv, DHSmalawi.csv, DHSmalawidst.csv, IGMEmalawi.csv, and RaMMPScalls.csv. These data files are not part of this repository but can be consolidated after running the Stata do-files RaMMPS_processing.do and IGME2024.do (also within this repository).

These routines require the following raw data files, available from these repositories:

MW_AnalyticSample.dta and MW_AllCallAttempts.dta from DataFirst, the RaMMPS project data repository, available at: https://doi.org/10.25828/M86Z-NF08

UN_IGME_2024.csv, the 2024 estimates from UN Inter-agency Group for Child Mortality Estimation, available at: https://childmortality.org/ 

MWIR7AFL.dta, MWBR7AFL.dta, and MWPR7AFL.dta, from the 2015–16 Malawi Demographic and Health Survey by the DHS Program, available at: https://dhsprogram.com

hh.sav, wm.sav, and bh.sav, from the 2019–20 Malawi Multiple Indicator Cluster Survey by UNICEF MICS, available at: https://mics.unicef.org

The m-file and do-files run automatically from top to bottom, but the user may need to adjust the file paths for reading the data and saving the outputs. The m-file requires some nested functions (also within this repository) to produce tables and some part of the analysis. 
# RaMMPS_U5M
