"""   
/   ********************************************************        
    * Metal Hydride ADC Rev 3	for cp2k		   *     
    ********************************************************  /  
"""


import pandas as pd
import numpy as np
import re

#user input
initfil = 'Gopt_pos-start.xyz'
endfil = 'Gopt_pos-end.xyz'
layers = 4
#

#Calculates average atomic displacement per layer.
def layripple():
	#grab data
	initconf = np.genfromtxt(initfil,delimiter=',') 
        endconf = np.genfromtxt(endfil,delimiter=',')
        
	
	geom0 = pd.DataFrame((initconf), columns=['u','v','w']) #list to dataframe
	geom1 = pd.DataFrame((endconf), columns=['u','v','w']) #list to dataframe


	#change dtype of columns
	geom0, geom1  = geom0.astype(float), geom1.astype(float)
	
	#convert pandas to numpy arrays
	geom0, geom1  = geom0.values, geom1.values
	
	#M and H seperation
	row, col = np.shape(geom0) #count rows, number will be used to seperate M and H

	#combine geom0 and geom1 into 1 array. needed to sort all values in same order
	both = np.hstack((geom0,geom1))

	#seperate into hydrogen containing geom0 amd geom1 and metal containing...
	hydrogen = both[0:(row/2),:]
	metal = both[(row/2):row,:]
	
	"""
	#argsort code example
	#ind = np.argsort(hydrogen[:,2]) #argsort is indirect sort
	#hydrogen = hydrogen[ind] #reconstruct array in sorted order
	"""
	
	#sort order by z component of geom0 (into layers)
	hydrogen = hydrogen[(np.argsort(hydrogen[:,2]))]
	metal = metal[(np.argsort(metal[:,2]))]

	#after sort seperate back into geom0 and geom1
	hydrogen0 = hydrogen[:,0:3]
	hydrogen1 = hydrogen[:,3:]
	metal0 = metal[:,0:3]
	metal1 = metal[:,3:]

	
	##start distance atoms move calculation

	hydrogen = np.absolute(hydrogen0 - hydrogen1)
	metal = np.absolute(metal0 - metal1)

	#hydrogen
	hydrogen_squared = hydrogen * hydrogen
	Hdist = np.sum(hydrogen_squared,axis=1)
	Hdist = np.sqrt((np.absolute(Hdist)))
	
	#metal
	metal_squared = metal * metal
	Mdist = np.sum(metal_squared,axis=1)
	Mdist = np.sqrt((np.absolute(Mdist)))
	
	##end distance atoms move calculation

	#detmine number of M or H per layer (typically 18)
	atoms_per_lay = len(Hdist)/layers
	

	#reshape arrays per layer to sum
	Hdist = Hdist.reshape(-1,atoms_per_lay)
	Mdist = Mdist.reshape(-1,atoms_per_lay)

	#determine avg. movement per lay
	Hdist = np.sum(Hdist, axis=1) / atoms_per_lay
	Mdist = np.sum(Mdist, axis=1) / atoms_per_lay
	
	print "metal moved (Ang) ="
	print Mdist
	print "hydrogen moved (Ang) ="
	print Hdist
	
	#MH pair movement	
	pairmov = Mdist + Hdist
	print "pair moved (Ang) ="
	print pairmov
