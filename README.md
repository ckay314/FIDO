# FIDO
GUI that takes in CME positional information and propagates a flux rope past a synthetic spacecraft and compares with observed data.  Details of the libraries needed to run the code, how to execute it, and a tutorial on how to begin understanding how the input parameters affect the profile can be found in FIDOmanual.pdf  
The ensemble version runs a set of cases with small, random variations in the CME latitude, longitude, and tilt.  This results in a pickle containing the resulting magnetic profiles (with rather coarse resolution).  The results can be displayed using insitusigplot, which compares with observations and determines the mean and one sigma range of the ensemble.

## FIDO-SIT
The FIDO code also contains the functionality to run FIDO-SIT, which includes the sheath.  The details of this are not fully explained in the manual but the sheath can be added by setting Add\_Sheath set to True in the input file, as well as the values for Sheath\_start, Sheath\_time, Compression, and Sheaeth\_v.  The Kp index can be added by setting Indices to True. 