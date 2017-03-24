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
        
    #store lines from file in data
    data = []

    for line_no in range (first_line,linecounter+1):
        grabline = linecache.getline(cifout, line_no)
        grabline = grabline.split() #comma seperate
        data.append(grabline)
    
    #convert to np
    data = np.array(data)
    return data[:,:-2]

#default angles, angles must be entered as createtable(['015','086']) etc.
def createtable(angles=['000','015','030','045','060',
                        '075','090','105','120','135',
                        '150','165','180']):

    file_count = len(angles) #number of files

    #1.Make A Blank Table
    header = np.array(['Angle', 'N_H1','N_H2','N_H3','Li_N','Li_H1','Li_H2','Li_H3','A'])
    header = header.reshape(1,9) #otherwise get shape (,9). numpy's error
    #make a blank 'nan' row for each angle.
    blankrows = np.full((file_count,9),np.nan)
    MasterTable = np.append(header,blankrows,axis=0)
    #!1
    

    for item in range(file_count):
        
        #2.Global Variables
        current_file = filegrabberCIF(angles[item]) #grabs cif
        current_row = item+1 #row number for current angle
        MasterTable[(current_row),0] = angles[item] #Col0 = current Angle
        #!2

        #create blank arrays to populate
        H_ammonia, N_ammonia = np.zeros(0), np.zeros(0)

        #3.Hydrogen in Ammonia
        #find Ammonia's Hydrogens
        for idxH in range(1,7): #(1,7) as only first 6H correspond to 2NH3, index of H
            hydrogen = np.where(current_file[:,0]=='H%d' % idxH) #find hydrogens 1 to 6 in current_file
            H_ammonia = np.append(H_ammonia, current_file[hydrogen])
            H_ammonia = np.reshape(H_ammonia,(-1,4))

        #determine equivalent hydrogens
        if np.shape(H_ammonia)[0] <3: #if less than 3H for ammonia in the CIF file, then some H must be equivalent
            if H_ammonia[1,0]=='H2': #then N--H2 and N--H3 are equiv.
                H_ammonia = np.vstack((H_ammonia,H_ammonia[1,:]))
            elif H_ammonia[1,0]=='H3': #then N--H1 and N--H2 are equiv.
                H_ammonia = np.vstack((H_ammonia[0,:],H_ammonia[1,:],H_ammonia[0,:]))
        #!3
        
        
        #4.Nitrogen in Ammonia
        for idxN in range(1,3):#(1,3) as 2N at most 
            nitrogen = np.where(current_file[:,0]=='N%d' % idxN)
            N_ammonia = np.append(N_ammonia, current_file[nitrogen])
            N_ammonia = np.reshape(N_ammonia,(-1,4))
        #!4
        
        #5.Find LithiumBelowAmmonia. *Assumption made that LithiumBelowAmmonia has the largest z-value of all the Lithiums.
        #find Li in current_file by locating rows containing substring 'Li'
        Li_substr = np.char.find(current_file[:,0],'Li') #returns 0 for rows containing 'Li'
        Lithium_rows = np.where(Li_substr ==0 ) #finds all 0 corresponding to 'Li'      
        LiFirstRow = np.min(Lithium_rows)
        LiLastRow = np.max(Lithium_rows) + 1
        
        #from current_file take lithium rows only
        Lithium = current_file[LiFirstRow:LiLastRow,:]
        Lithium = np.reshape(Lithium,(-1,4))
        
        #CIF FILE --> to find Lithium atom furthest from surface need to also check inverted positions because its a cif file
        InvLithium = 1 - (Lithium[:,1:].astype(np.float)) #inverting XYZ coordinates and not column idx
        InvLithium = np.hstack(((Lithium[:,0].reshape(-1,1)),
                                InvLithium)) #stack index column to inverted XYZ columns

        #Stack ALL Lithiums
        LithiumStack = np.vstack((Lithium,InvLithium))
        #find Lithium Atom with largest z value
        LithiumZStack = LithiumStack[:,3].astype(float) #z col. only
        #find row with largest z-value from LithiumZStack
        LithiumBelowAmmonia = LithiumStack[np.where(LithiumZStack == np.max(LithiumZStack)),:]       
        LithiumBelowAmmonia = np.reshape(LithiumBelowAmmonia, (-1,4))

        """
        #Manually overide LithiumBelowAmmonia Position for 5 Ang calc
        findLirow = np.where(Lithium =='Li16')[0]
        LithiumBelowAmmonia = Lithium[findLirow,:]
        """
        #!5

        #6.Convert rows of hydrogen, nitrogen and lithium selected from current_file into a usable format
        #remove index and convert to float
        N_ammonia = N_ammonia[:,1:].astype(float) 
        H_ammonia = H_ammonia[:,1:].astype(float)     
        LithiumBelowAmmonia = LithiumBelowAmmonia[:,1:].astype(float)
 
        #convert frac_coord to absolute positions
        H_ammonia[:,:2] *= ab
        H_ammonia[:,2] *= c
        N_ammonia[:,:2] *= ab
        N_ammonia[:,2] *= c
        LithiumBelowAmmonia[:,:2] *=ab
        LithiumBelowAmmonia[:,2] *=c
        #!6
        
        #7.Number of rows
        H_rows = np.shape(H_ammonia)[0]
        N_rows = np.shape(N_ammonia)[0]
        Li_rows = np.shape(LithiumBelowAmmonia)[0]
        #!7

        #8.Internal SymmetryCheck
        #if too many rows then symmetry has failed
        if H_rows >3 or N_rows >1 or Li_rows > 1:
            return 'InvSymm Failed'
        #!8
        
        #9.Calculate Bond Lengths
        #bond lengths involving hydrogen, iterating over 3H
        for i in range(H_rows):
            #nitrogen-hydrogen bond length
            AmmoniaBL = np.sqrt(np.sum( ((H_ammonia[i] - N_ammonia)**2), axis=1)) #i here is row
            #lithium-hydrogen bond length
            Li__H3N_BL = np.sqrt(np.sum( ((H_ammonia[i] - LithiumBelowAmmonia)**2), axis=1))
            
            #add to table
            MasterTable[current_row,(i+1)] = float(AmmoniaBL)
            MasterTable[current_row,(i+5)] = float(Li__H3N_BL) #i+5 is column
    
        #lithium-nitrogen bond length
        LiN_BL = np.sqrt(np.sum( ((N_ammonia - LithiumBelowAmmonia)**2)))
        MasterTable[current_row, 4] = float(LiN_BL)
        #!9
    
    return MasterTable[1:,:]
