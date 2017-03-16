import pandas as pd
import numpy as np
import linecache

def filegrabber(angle):
    castepf = 'GOPT_L-BFGS_LiH33_11LAY_750eV_kpn551_wFixedBulk_44AngBox_SymmGen_w2NH3at_4AngfromSurf_%sdeg.castep' % angle

    #locates line of geom1 in castepfile
    with open(castepf, 'r') as myfile:
        for (i, line) in enumerate(myfile):
            if 'Total number of ions in cell =' in line:
                #returns last word in the line
                no_of_ions = line.split()[-1]
            #locates lines in castep file of geo
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
    #has to be at this stage to calc distance between atoms
    #can't calculate dist between atoms in frac. coords
    ab = 12.0251939
    c  = 44.0420257 
    data[:,:2] *=ab
    data[:,2] *=c
    return data

#default angles, angles must be entered as createtable(['015','086']) etc.
def createtable(angles=['000','015','030','045','060',
                        '075','090','105','120','135',
                        '150','165','180']
               ):
   
    #loops over filegrabber() for every file
    for item in range(len(angles)):
        #setup blank table to populate
        table = np.zeros((1,8))
        
        current_file = filegrabber(angles[item])

        #pass current_file to NH bond calculator
        AmmoniaBonds(current_file)
        
        #pass to NLi length
        NLi = DistNLi(current_file)
        print NLi
               

        #output table as txtfile

def AmmoniaBonds(current_file):
    #grab xyz of N from current_file
    for N_rowno in range(1,3):
        N_atoms = current_file[-(N_rowno)]
        
        #poor solution, update in later version
        if N_rowno == 1:
            Counter = 3
        else:
            Counter = 0
        
        #grab 3H atoms associated with N
        H_atoms = current_file[0+Counter:3+Counter]
        
        #calculate dist between atoms
        d = np.sqrt(np.sum( ((H_atoms - N_atoms)**2)  ,axis=1))
        print d

def DistNLi(current_file):
    #Specific xyz coordinates
    N_Bottom = current_file[-1]
    Li_Bottom = current_file[219]    
    N_Top = current_file[-2]
    Li_Top =  current_file[227]

    #calculate distance
    NLi_Top = np.sqrt(np.sum( ((N_Top - Li_Top)**2)))
    NLi_Bottom = np.sqrt(np.sum( ((N_Bottom - Li_Bottom)**2)))
    
    return [NLi_Top, NLi_Bottom]
