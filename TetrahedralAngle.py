import numpy as np
import pandas as pd

outfile = "out"
infile = "in.xyz"

def timestep(angle):
    #print "time steps run t=0 to t=odd number as gives an even number of rows"
    #print "ANY error related array dimensions, delete last timestep"
    row_count = len(angle) / 2
    row_count = row_count * 0.001 #fs to ps

    time = np.arange(0,row_count,0.0005)
    return np.column_stack((time,angle))

def Get_Coordinates():
	row_count = 510 #number of atoms + header=
	PointsofTetrahedron = [504,505,506,507] #[500,501,502,503] ammonia atomic positions 

	data = pd.read_csv(infile, skiprows=2, header=None,\
			   error_bad_lines=False, names=['id', 'x','y','z'], sep="\s+")

	CoordArray = np.array([])

	for row_id in PointsofTetrahedron:
		Atom_coord = data.iloc[row_id::row_count,1:].astype('float')
        	Atom_coord = Atom_coord.values #to numpy
		CoordArray = np.append(CoordArray, Atom_coord) #append all coordinates to a single array
	
	#split into sub arrays
	A, B, C, D = np.split(CoordArray, 4,axis=0)
	#have to reshape now, can't be done earlier
	A,B,C,D = A.reshape(-1,3), B.reshape(-1,3),\
		  C.reshape(-1,3), D.reshape(-1,3)

	return A,B,C,D

def CentroidofaTetrahedron(rebin):
	#!Cal_CentroidofaTetrahedron
	A, B, C, D = Get_Coordinates()
	Centroid_x = (A[:,0] + B[:,0] + C[:,0] + D[:,0]) / 4.0
	Centroid_y = (A[:,1] + B[:,1] + C[:,1] + D[:,1]) / 4.0
	Centroid_z = (A[:,2] + B[:,2] + C[:,2] + D[:,2]) / 4.0
	PointC = np.stack((Centroid_x, Centroid_y, Centroid_z), axis=1) #centre point
	#!Cal_CentroidofaTetrahedron

	#!ABCangle
	#Use np.copy as otherwise PointB = PointA.
	#Otherwise any modification to PointA would modify PointB, as are both coupled to variable D.	
	#https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.copy.html
	PointA = np.copy(D)
	PointB = np.copy(D)
	PointA[:,2] += 5

	vectorAB = PointB - PointA
	vectorBC = PointC - PointB

	dotproduct = (vectorAB[:,0] * vectorBC[:,0]) +\
                     (vectorAB[:,1] * vectorBC[:,1]) +\
                     (vectorAB[:,2] * vectorBC[:,2])
        
	pythagAB = np.sqrt((vectorAB[:,0]**2) +\
                           (vectorAB[:,1]**2) +\
                           (vectorAB[:,2]**2))
                
        pythagBC = np.sqrt((vectorBC[:,0]**2) +\
                           (vectorBC[:,1]**2) +\
                           (vectorBC[:,2]**2))

        angle = np.arccos((-1*dotproduct)/(pythagAB * pythagBC))
	#angle = np.degrees(angle) convert to degrees
	#!ABCangle

	timestep_angle = timestep(angle) #add time column
	x, y = timestep_angle[:,0], timestep_angle[:,1] #seperate into x and y
	
	#!Rebin
	#Modified from https://stackoverflow.com/questions/30379311/fast-way-to-take-average-of-every-n-rows-in-a-npy-array
	x, y = x.reshape(1,-1), y.reshape(1,-1)
	x = x.transpose().reshape(-1,rebin).mean(1).reshape(1,-1).transpose()
	y = y.transpose().reshape(-1,rebin).mean(1).reshape(1,-1).transpose()
	#!rebin

	#return np.savetxt(outfile, timestep_angle)
	return x, y

def Plot():
	import matplotlib as mpl
	import matplotlib.pyplot as plt

	x, y = CentroidofaTetrahedron(rebin=1)
	x_rebin,y_rebin = CentroidofaTetrahedron(rebin=20) #rebin has to be an even number

	mpl.rcParams['font.size']=20
	fig, ax1 = plt.subplots()
	fig.set_size_inches(40, 12)

	#!ax1    
	ax1.set_xlabel('time (ps)')	    
	ax1.set_ylabel('Orientation (rad)')

	ax1.set_ylim([0,3])
	ax1.set_xlim([0,8])

	ax1.set_yticks([0.,  0.16666666667*np.pi,\
			0.33333333333*np.pi, 0.5*np.pi,\
			0.66666666667*np.pi, 0.83333333333*np.pi,\
			np.pi])

	ax1.set_yticklabels([u'up \u2192 0',  r"$\frac{\pi}{6}$",\
			     r"$\frac{\pi}{3}$", r"$\frac{\pi}{2}$",\
			     r"$\frac{2\pi}{3}$",  r"$\frac{5\pi}{6}$",\
			     u"down \u2192 $\pi$"])

	ax1.tick_params(axis='both', direction="in",length=10,top=True,bottom=True,left=True)

	ax1.plot(x,y, 'b', lw=4, alpha=0.3, label="Original")
	ax1.plot(x_rebin,y_rebin, 'b', lw=2.5, label="Rebin = 20")
	#!ax1    

	#!ax2_data
	A, B, C, D = Get_Coordinates() #Suboptimal re-calling function but still fast.
	bulk = 18.040304 #relaxed bulk value
	NitrogenHeightAbove = D[:,2] - bulk #z_component of Nitrogen Above Surface
	NitrogenHeightAbove = timestep(NitrogenHeightAbove)
	x_2, y_2 = NitrogenHeightAbove[:,0], NitrogenHeightAbove[:,1]
	#!ax2_data

	#ax2
	ax2 = ax1.twinx()
        ax2.tick_params(direction="in", length=10)
	ax2.set_ylabel(r'height from surface ($\AA$)')
	ax2.set_ylim([0,5])
	ax1.plot(np.nan,'b', lw=2, alpha=0.2, label = 'Approximate height')  # to append ax2 label to ax1 legend
	ax1.legend(loc=0) #Need location for legend to show. Silly.
	ax2.plot(x_2,y_2, 'b', lw=2, alpha=0.2)
	#!ax2

	#!GreyBoxes
	"""	
	#1NH3 case
	ax1.axvspan(0,0.75, color='k', alpha=0.1)
	ax1.axvspan(2.3565,3.865, color='k', alpha=0.1)
	"""
	
	"""
	#2NH3 case
	#atoms [500,501,502,503] 
	ax1.axvspan(0,0.989, color='k', alpha=0.1)
	ax1.axvspan(17.5,20, color='k',alpha=0.1)
	"""
	#!GrayBoxes

	fig.show()
	fig.savefig(outfile, dpi=300)
