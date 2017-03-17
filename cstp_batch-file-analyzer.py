#note: Li--N Code is not robust enough and needs commenting

import pandas as pd
import numpy as np
import linecache

#unit cell
ab = 12.0251939
c  = 44.0420257

#poor solution which needs improvement
def filegrabberCIF(angle):
    #use castep file not cif
    cifout = 'GOPT_L-BFGS_LiH33_11LAY_750eV_kpn551_wFixedBulk_44AngBox_SymmGen_w2NH3at_4AngfromSurf_%sdeg-out.cif' % angle
    
    linecounter = 0
   
    with open(cifout, 'r') as myfile:
        for (line_no, line) in enumerate(myfile):

            #finds line where coordinates start
            if '_atom_site_occupancy' in line:
                first_line = line_no + 2
            linecounter += 1
        
    #store lines in rows
    rows = []

    for line_no in range (first_line,linecounter+1):
        grabline = linecache.getline(cifout, line_no)
        grabline = grabline.split() #comma seperate
        rows.append(grabline)
    
    #convert to np
    rows = np.array(rows)
    return rows[:,:-2]

#default angles, angles must be entered as createtable(['015','086']) etc.
def createtable(angles=['000','015','030','045','060',
                        '075','090','105','120','135',
                        '150','165','180']
               ):
   
    #loops over filegrabber() for every file
    for item in range(len(angles)):
                
        current_file = filegrabberCIF(angles[item])
        
        #make blank np arrays
        AmmoniasH = np.zeros(0)
        AmmoniasN = np.zeros(0)

        #finds H in Ammonia
        for H_idx in range(1,7):
            Hrow = np.where(current_file[:,0]=='H%d' % H_idx)
            AmmoniasH = np.append(AmmoniasH, current_file[Hrow])

        for N_idx in range(1,3):
            Nrow = np.where(current_file[:,0]=='N%d' % N_idx)
            AmmoniasN = np.append(AmmoniasN, current_file[Nrow])
        
        #Lithium atom below N
        LiBelowAmmonia = current_file[70]
        print LiBelowAmmonia

        #reshape
        AmmoniasH = np.reshape(AmmoniasH,(-1,4)) 
        AmmoniasN = np.reshape(AmmoniasN,(-1,4))
        LiBelowAmmonia = np.reshape(LiBelowAmmonia,(-1,4))

        #remove indx and convert to float
        AmmoniasH = AmmoniasH[:,1:].astype(float)
        AmmoniasN = AmmoniasN[:,1:].astype(float) 
        LiBelowAmmonia = LiBelowAmmonia[:,1:].astype(float)

        #LiCode
        #check Li is the correct one!
        print LiBelowAmmonia[0,2]
        print AmmoniasN[0,2]
        #if surf Li and N are greater than ~13Ang apart 
        if abs(LiBelowAmmonia[0,2]-AmmoniasN[0,2]) > 0.3:
            LiBelowAmmonia = 1 - LiBelowAmmonia
        print LiBelowAmmonia

        
        #convert fraccoord to absolute positions
        AmmoniasH[:,:2] *= ab
        AmmoniasN[:,:2] *= ab
        LiBelowAmmonia[:,:2] *=ab

        AmmoniasH[:,2] *= c
        AmmoniasN[:,2] *= c
        LiBelowAmmonia[:,2] *=c


        #calc num_rows
        num_Hrows = np.shape(AmmoniasH)[0]
        num_Nrows = np.shape(AmmoniasN)[0]
        num_Lirows = np.shape(LiBelowAmmonia)[0]

        #SYMMETRYCHECKER
        if num_Hrows >3 or num_Nrows >1 or num_Lirows > 1:
            return 'InvSymm Failed'
        #SYMMETRYCHECKER
        
        #elements coordinates to functions, comment out as appropriate
        AmmoniaBL(AmmoniasN,AmmoniasH,num_Hrows)    
        LiN_BL(AmmoniasN, LiBelowAmmonia) 
        #LiH_BL(AmmoniasH, LiBelowAmmonia)

#Ammonia Bond Lengths
def AmmoniaBL(AmmoniasN,AmmoniasH,num_Hrows):
        
    #grab coordinates of each ammonia in turn
    for rows in range(num_Hrows):        
        #d is dist btween atoms
        d = np.sqrt(np.sum( ((AmmoniasH[rows] - AmmoniasN)**2)  ,axis=1))
        print 'N--H%d = ' % (rows +1)
        print d

def LiN_BL(AmmoniasN, LiBelowAmmonia):
 
    LiN_BondLength = np.sqrt(np.sum( ((AmmoniasN - LiBelowAmmonia)**2)))
    
    print 'Li--N='
    print LiN_BondLength
