import pandas as pd
import numpy as np
import linecache

def filegrab():
    castepf = 'GOPT_L-BFGS_LiH33_11LAYERS_wFixedBulk_750eV_w2NH3at_2AngfromSurf_000deg.castep'

    #locates line of geom1 in castepfile
    with open(castepf, 'r') as myfile:
        for (i, line) in enumerate(myfile):
            if 'Total number of ions in cell =' in line:
                #returns last word in the line
                no_of_ions = line.split()[-1]
            #locates lines in castep file of geom1
            if 'Final Configuration' in line:
                lineofgeom1 = i + 12

    #https://www.safaribooksonline.com/library/view/python-cookbook-2nd/0596007973/ch02s05.html
    datalist=[]

    for m in range (lineofgeom1,lineofgeom1 + int(no_of_ions)):
        theline = linecache.getline(castepf, m)
        theline = theline.split() #comma seperate
        theline.pop(0) ,theline.pop(-1) #removes 'x' from beginning and end
        datalist.append(theline)
    
    #convert to np and remove first two cols
    data = np.array(datalist)
    data = data[:,2:]
    data = data.astype(float)

    #convert to absolute positions
    ab = 12.0251939
    c  = 44.0420257 
    data[:,:2] *=ab
    data[:,2] *=c

    return data

def gentable():
    cfile = filegrab()
    
    
    #setup blank table to populate
    table = np.zeros((1,8))

    #Calc NH bond lengths
    for N_rowno in range(1,3):
        NitrogenAtm = cfile[-(N_rowno)]
        
        #simple fix
        if N_rowno == 1:
            n = 3
        else:
            n = 0

        hydrogenAtm = cfile[0+n:3+n]
        d = (hydrogenAtm - NitrogenAtm)**2
        d = np.sqrt(np.sum(d,axis=1))
        print d
