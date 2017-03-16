#have to choose index (uncomment the one you want torun) and manually fill in ions_to_fix

import numpy as np
import pandas as pd

#generate ionic constraints
def genconst(element):  #enter 'Li' not Li    
    #specify positions

    ions_to_fix =[151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168]
    #ions_to_fix =[181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198]
    #cartesian coordinates to constrain
    constr_coord  = np.array([[[1,0,0],
                        [0,1,0],
                        [0,0,1]]])

    constr_coord = np.repeat(constr_coord,(len(ions_to_fix)), axis=0)
    constr_coord = constr_coord.reshape(-1,3)
    
    ions_to_fix = np.repeat(ions_to_fix,3)
    ions_to_fix = ions_to_fix.reshape(1,-1)
    ions_to_fix = ions_to_fix.T
    
    
    ionic_constr = np.concatenate((ions_to_fix,constr_coord),
                                  axis=1)
   
    #making an index
    #obj=0 meaning insert idx to first column in ionic_constr
    idx_length = len(ionic_constr)

    #choose one idx
    idx = [i for i in range(1, (idx_length+1))]
    #idx = [i for i in range((idx_length+1), ((2*idx_length)+1) )]

    ionic_constr = np.insert(ionic_constr,0, idx, axis=1)

    #element column
    el = []
    for m in range(1, len(ionic_constr)+1):
        el.append(element)
    el = pd.Series(el,name='Element')


   
    #convert numpy array to Dataframe so a column of strings can be added
    ionic_constr =  pd.DataFrame(ionic_constr, columns=['idx','ionPos','u','w','v'])
    #concat elements (Series) to ionic_constr (Dataframe)
    ionic_constr = pd.concat([el,ionic_constr],axis=1)
    
    #reorder dataframe columns
    ionic_constr = ionic_constr[['idx','Element','ionPos','u','w','v']]
    
    ionic_constr.to_csv('/home/mimas/Desktop/ionic_constr.cell',sep='\t',header=False,index=False)
