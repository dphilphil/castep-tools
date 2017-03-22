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

    #make table
    header = np.array(['Angle', 'N_H1','N_H2','N_H3','Li_N','Li_H1','Li_H2'])
    header = header.reshape(1,7) #otherwise get shape (,7). numpy's error
    
    #generate table of correct length    
    rowstoadd = len(angles)
    #make blankrows
    blankrows = np.full((rowstoadd,7),np.nan)
    table = np.append(header,blankrows,axis=0)

    #loops over filegrabber() for every file
    for item in range(len(angles)):
        g_row = item+1 #global row
        table[(g_row),0] = angles[item]
        current_file = filegrabberCIF(angles[item])
        
        #make blank np arrays
        AmmoniasH, AmmoniasN = np.zeros(0), np.zeros(0)

        #finds H in Ammonia
        for H_idx in range(1,7):
            Hrow = np.where(current_file[:,0]=='H%d' % H_idx)
            AmmoniasH = np.append(AmmoniasH, current_file[Hrow])
        
        AmmoniasH = np.reshape(AmmoniasH,(-1,4))

        if np.shape(AmmoniasH)[0] <3: #if less than 3 H atoms for NH3 in cif cell
            if AmmoniasH[1,0]=='H2':
                AmmoniasH = np.vstack((AmmoniasH,AmmoniasH[1,:]))
            elif AmmoniasH[1,0]=='H3':
                AmmoniasH = np.vstack((AmmoniasH[0,:],
                                       AmmoniasH[1,:],
                                       AmmoniasH[0,:],
                                      ))
        #print AmmoniasH
        
        
        for N_idx in range(1,3):
            Nrow = np.where(current_file[:,0]=='N%d' % N_idx)
            AmmoniasN = np.append(AmmoniasN, current_file[Nrow])
        
        #!!!!!!!!!!!!!!!!!!!!LI!!!!!!!!!!!!!!
        
        #grab all Li atoms
        #find Li with largest z
        #locate rows containing substring 'Li'
        li_substr = np.char.find(current_file[:,0],'Li')
        #find Li index
        Li_idx = np.where(li_substr ==0 )
        first_row = np.min(Li_idx)
        last_row = (np.max(Li_idx)) +1
        #grab all Li
        All_Li = current_file[first_row:last_row,:]
      
        All_Li = np.reshape(All_Li,(-1,4)) #code needed!

        #np.append(All_Li,All_Li)
        #b = All_Li[:,1:].astype(np.float)

        #creating a table of inverted Li positions
        invertedLixyz = 1 - (All_Li[:,1:].astype(np.float))
        LiCol1 = All_Li[:,0].reshape(-1,1) #have to reshape!
        
        invertedLi = np.hstack((LiCol1,invertedLixyz))
        
        #stack Li and Li Inverted
        All_Li = np.vstack((All_Li,invertedLi))

        """code for finding Li with largest Z component"""
        ZLi = All_Li[:,3].astype(float) 
        #find idx of Li with greatest z value, i.e LibelowN
        findLirow = np.where(ZLi == np.max(ZLi))
        LiBelowAmmonia = All_Li[findLirow,:]
        LiBelowAmmonia = np.reshape(LiBelowAmmonia, (-1,4))

        """
        #manually overide LiAtom Position for 4Ang and 5 Ang
        findLirow = np.where(All_Li =='Li16')[0]
        LiBelowAmmonia = All_Li[findLirow,:]
        """

        #reshape
        AmmoniasN = np.reshape(AmmoniasN,(-1,4))

        #remove indx and convert to float
        AmmoniasH = AmmoniasH[:,1:].astype(float)
        AmmoniasN = AmmoniasN[:,1:].astype(float) 
        LiBelowAmmonia = LiBelowAmmonia[:,1:].astype(float)

        #ensures that Li and N atoms are on the sameside
        #Li determined to be on the otherside to N if Li--N distance is > ~13Ang 
        if abs(LiBelowAmmonia[0,2]-AmmoniasN[0,2]) > 0.3:
            LiBelowAmmonia = 1 - LiBelowAmmonia
        #print LiBelowAmmonia

        
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
        
        #Ammonia Bond Length
        #grab coordinates of each ammonia in turn
        for row in range(num_Hrows):
            AmmoniaBL = np.sqrt(np.sum( ((AmmoniasH[row] - AmmoniasN)**2)  ,axis=1))
            #print 'N--H%d = ' % (row +1)
            #print AmmoniaBL
            table[g_row,(row+1)] = float(AmmoniaBL)
        #Ammonia Bond Length
    
        #Lithium nitrogen Bond Length
        LiN_BL = np.sqrt(np.sum( ((AmmoniasN - LiBelowAmmonia)**2)))
        #print 'Li--N=' + str(LiN_BL)
        table[g_row, 4] = float(LiN_BL)
    
    print table
    

def Li_H3N_BL(AmmoniasH, LiBelowAmmonia,num_Hrows):
    for rows in range(num_Hrows):
        d = np.sqrt(np.sum( ((AmmoniasH[rows] - LiBelowAmmonia)**2)  ,axis=1))
        print 'Li--H%d = ' % (rows +1)
        print d
