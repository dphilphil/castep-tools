"""   
/   ******************************************************** 
    *                                                      *        
    * Metal Hydride ADC Rev 2				   *   
    * Author:  Phillip Marks                               *   
    * Date:    20-09-2017                                  *     
    * Source:  https://github.com/dphilphil/castepy/       *   
    *                                                      *   
    ********************************************************  /  
"""


import pandas as pd
import numpy as np
import re

#user input
initfil = 'Gopt_KH333_7LaysRep13_fromGoptKH111_18AngV_2BottomLayFixed_0.01gft.cell'
endfil = 'Gopt_KH333_7LaysRep13_fromGoptKH111_18AngV_2BottomLayFixed_0.01gft-out.cell'
ab = 17.07439206641234
c = 35.07439206641234
layers = 7
#

def initconf():
	data = [] #blank list

	with open(initfil,'r') as fileinput:
		for line in fileinput:
			line = line.lower() #make lowercase
			if line.strip() == '%block positions_frac':
				break
		for line in fileinput:
			line = line.lower()
			if line.strip() == '%endblock positions_frac':
				break

			line = line.rstrip('\n') #remove \n from line
			line = re.split(r'\t+', line) #split line based on \t
			
        		data.append(line)

		return data


def endconf():
	data = [] #blank list

	with open(endfil,'r') as fileinput:
		for line in fileinput:
			line = line.lower() #make lowercase
			if line.strip() == '%block positions_frac':
				break
		for line in fileinput:
			line = line.lower()
			if line.strip() == '%endblock positions_frac':
				break
			
			line = line.rstrip("\n") #remove \n from line
			line = line.lstrip(' ') #removes leading whitespace from start of lines
			line = re.split("\s\s+" , line) #split based on whitespace
			
        		data.append(line)

		return data

#only works without NH3
#Calculates average atomic displacement per layer.
def layripple():

	geom0 = pd.DataFrame((initconf()), columns=['element','u','v','w']) #list to dataframe
	geom1 = pd.DataFrame((endconf()), columns=['element','u','v','w']) #list to dataframe
	
	#delete element column
	geom0, geom1 = geom0.drop('element', 1), geom1.drop('element', 1)

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
	hydrogen[:,0:2] = (hydrogen[:,0:2]) * ab	#multiply cols 0 and 1 by ab 
	hydrogen[:,2] = (hydrogen[:,2]) * c
	hydrogen_squared = hydrogen * hydrogen
	Hdist = np.sum(hydrogen_squared,axis=1)
	Hdist = np.sqrt((np.absolute(Hdist)))
	
	#metal
	metal[:,0:2] = (metal[:,0:2]) * ab
	metal[:,2] = (metal[:,2]) * c
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
