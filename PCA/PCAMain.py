# Author: Zhinuo J. Wang 5 July 2016
#####################################################################################
# Use principal component analysis (PCA) to extract geometric modes
# and modal coefficients for a population of left ventricular (LV) models 
# segmented from patient MRI. These modes and modal coefficients will then be used
# in a framework for estimating the reference geometry (i.e. estimating the
# modal coefficients which combine to re-create the LV geometry) as well as the
# passive stiffness parameter simultaneously. 

# Input(s):
# 	- 28 patient LV models in 16 tricubic Hermite elements. 
# Output(s):
# 	- All geometric modes which describe variation in geometry within the 
#	  popluation.
#	- Corresponding eigenvalues/modal coefficients. 
# 	- Small number of modes capable of describing 99% of the geometric variation. 
#####################################################################################



