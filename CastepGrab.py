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

    return data

def gentable():
    dj = filegrab()
    #convert to float
    dj = dj.astype(float)
    #setup blank table to populate
    table = np.zeros((1,8))
    
    for N_rowno in range(-2,0):
        NitrogenAtm = dj[N_rowno]
        print NitrogenAtm
        
        #calculate NH bond lengths
        for rowno in range(6):
            #grab hydrogen atom coords
            hydrogenAtm = dj[rowno]
            print hydrogenAtm
            #(x0,y0,z0), (x1,y1,z1), d = sqrt( (x1-x0)^2 + (y1-y0)^2 + (z1-z0)^2 )
            print hydrogenAtm[0]
            #d = np.sqrt ( (hydrogenAtm[0]
