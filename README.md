# FIDO
GUI that takes in CME positional information and propagates a flux rope past a synthetic spacecraft and compares with observed data.  Old versions hardcoded to pull in information from a pickle about one of 45 CMEs from Kay et al. 2017 but new can read in details from a formatted text file.  
The ensemble version runs a set of cases with small, random variations in the CME latitude, longitude, and tilt.  This results in a pickle containing the resulting magnetic profiles (with rather coarse resolution).  The results can be displayed using insitusigplot, which compares with observations and determines the mean and one sigma range of the ensemble.