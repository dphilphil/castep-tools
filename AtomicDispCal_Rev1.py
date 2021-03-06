"""   
/   ******************************************************** 
    *                                                      *        
    * Metal Hydride ADC Rev 1				                    *   
    * Author:  Phillip Marks                               *   
    * Date:    25-04-2017                                  *     
    * Source:  https://github.com/dphilphil/castepy/       *   
    *                                                      *   
    ********************************************************  /  
"""

import pandas as pd
import numpy as np
import linecache

#fname use %s
fname = 'GOPT_L-BFGS_LiH33_%sLAYrepresenting11LAY_wFixed2BottomLays_750eV_kpn551_28AngBox_dipolecorr_NH3at2.6Ang_%sdeg.castep'
#fname = 'GOPT_L-BFGS_LiH33_%sLAYrepresenting11LAY_wFixed2BottomLays_fromGOPTLiHv2_600eV_kpn551_NoSymmGen_12AngV_0.02ft.castep'
#fname = 'GOPT_L-BFGS_LiH33_6LAYrepresenting11LAY_wFixedBottom_fromGOPTLiHv2_600eV_kpn551_NoSymmGen_12AngV_0.02evang.castep'
#fname = 'GOPT_L-BFGS_LiH33_11LAYERS_wFixedBulk_fromGOPTLiHv2_600eV_kpn551_44AngV_SYMMGEN_0.02ft.castep'

def infile(layerno, Angle):
    castepf = fname % (layerno,Angle)

    #locates line of geom0 in castepfile
    with open(castepf, 'r') as myfile:
        for (i, line) in enumerate(myfile):
            if 'Total number of ions in cell =' in line:
                #returns last word in the line
                no_of_ions = line.split()[-1]
            #locates lines in castep file of geom0
            if 'Unit Cell' in line:
                lineofgeom0 = i + 26
    
    data_list=[]

    for m in range (lineofgeom0,lineofgeom0 + int(no_of_ions)):
        theline = linecache.getline(castepf, m)
        theline = theline.split() #comma seperate
        theline.pop(0) ,theline.pop(-1) #removes 'x' from beginning and end
        data_list.append(theline)

    return data_list

def outfile(layerno, Angle):
    castepf = fname % (layerno,Angle)

    #locates line of geom1 in castepfile
    with open(castepf, 'r') as myfile:
        for (i, line) in enumerate(myfile):
            if 'Total number of ions in cell =' in line:
                #returns last word in the line
                no_of_ions = line.split()[-1]
            #locates lines in castep file of geom1
            if 'Final Configuration' in line:
                lineofgeom1 = i + 12

    data_list=[]

    for m in range (lineofgeom1,lineofgeom1 + int(no_of_ions)):
        theline = linecache.getline(castepf, m)
        theline = theline.split() #comma seperate
        theline.pop(0) ,theline.pop(-1) #removes 'x' from beginning and end
        data_list.append(theline)
        
    return data_list

def halfcelldisp(angles=['000','060','090']):

    layerno = 6 #as halfcell of 11L LiH
    ab = 12.0251939
    c = 28.0210306

    #run code for each angle
    for AngleIdx in range (len(angles)):
        #individaully select angle
        Angle = angles[AngleIdx]
        
        #grab start and end configs.
        geom0 = np.array(infile(layerno,Angle))
        geom1 = np.array(outfile(layerno,Angle))
        #identifier contains element and indexes
        id = geom0[:,:2]
        #combine geom0 and geom geom1 into a single table only for sorting
        botharr = np.column_stack((id,geom0[:,2:],geom1[:,2:]))
        botharr = botharr[botharr[:,4].argsort()] #crudely sort by the z column of geom0
        #having sorted the table seperate components 
        sortedid = botharr[:,:2]
        sortedgeom0 = botharr[:,2:5].astype(float)
        sortedgeom1 = botharr[:,5:8].astype(float)
        #have to sort order then convert factional coordinates to absolute positions. argsort bug.    
        sortedgeom0[:,:2] *=ab
        sortedgeom0[:,2] *= c
        sortedgeom1[:,:2] *= ab
        sortedgeom1[:,2] *=c
        
        AtomicDisp = np.sqrt(np.sum(((sortedgeom1-sortedgeom0)**2),axis=1 )) #axis=1 as summing across each row 
        
        #for convenience only
        AtomicDisp = np.round(AtomicDisp,5)
        
        #DIRECTION OF VECTOR
        z_sign = np.sign(sortedgeom1[:,2] - sortedgeom0[:,2]) #sign of delta z used to give indication of the direction atoms move in
        AtomicDisp = np.column_stack((sortedid,(AtomicDisp*z_sign))) #(AtomicDisp*z_sign) - sign only used to denote direction in z
        
        #make seperate numpy tables for Elements 
        LiRows = np.where(AtomicDisp[:,0]=='Li') #finds all rows where first column contains 'Li'
        LithiumArr = AtomicDisp[LiRows,:]  #grab rows using index
        #combining two lines above
        HydrogenArr = AtomicDisp[(np.where(AtomicDisp[:,0]=='H')),:]      
        
        #combined elemental arrays only for the purpose of saving a table
        #CombinedArrayZero as contains rows where atoms hadn't moved i.e. 0Ang movement
        CombinedArr0 = np.hstack((LithiumArr,HydrogenArr))
        CombinedArr0 = np.reshape(CombinedArr0,(-1,3))    #reshape sois a single arr rather then each line being an array
        CombinedArr0 = np.insert(CombinedArr0,0,['E','id', (Angle+'deg') ],axis=0) #adds header
       
        #deleting rows where atoms hadn't moved as were fixed
        CombinedArr = CombinedArr0[(np.where(CombinedArr0[:,2]!='0.0')),:]      
        CombinedArr = np.reshape(CombinedArr, (-1,3)) 


        #generate Table containing data for all angles
        if AngleIdx==0:
            #for first angle in list (Angle Index = 0) genrate array
            Tab = CombinedArr
        else:
            #for subsequent angles in list, append data to table
            Tab = np.column_stack((Tab,CombinedArr[:,2])) #(AtomicDisp*z_sign) - sign only used to denote direction in z
        
    np.savetxt(fname +'_TABLE',Tab,
               fmt='%s',delimiter='\t'
              )

#Calculates average LiH atomic displacement per layer.
def LiH_rippling(layerno):
    
    geom0 = pd.DataFrame((infile(layerno)),
            columns=['Element','ionNumber','u','v','w'])

    #change dtype of columns
    geom0[['ionNumber', 'u', 'v','w']] = geom0[['ionNumber','u','v','w']].astype(float)

    """ if the order of table is not H Li H Li etc.
    then you have to round to fewer d.p. """
    
    #must round w as starting positions not all the same
    geom0['w'] = geom0['w'].round(3)
    geom0['u'] = geom0['u'].round(4)
    geom0['v'] = geom0['v'].round(4)
    
    #order geom0 based on w then u then v
    geom0 = geom0.sort_values(['w','u','v'],
                              ascending=[True,True,True])

    #find index order for sorted geom0
    geom0index = geom0.index.tolist()
    
    #geom1 is compared to geom0
    geom1 = pd.DataFrame((outfile(layerno)),
            columns=['Element','ionNumber','u','v','w'])

    geom1[['ionNumber', 'u', 'v','w']] = geom1[['ionNumber','u','v','w']].astype(float)
    
    #redorder geom1 using geom0index
    geom1 = geom1.reindex(geom0index)
   
    #!!have to get geom0 again w/o rounding!!
    geom0 = pd.DataFrame((infile(layerno)),
            columns=['Element','ionNumber','u','v','w'])


    geom0[['ionNumber', 'u', 'v','w']] = geom0[['ionNumber','u','v','w']].astype(float)
    
    geom0 = geom0.reindex(geom0index)
    #end of !!getting geom0 again!!

    #!!!!!!!!this is where code for Li or H would go!!!
    #geom0 = geom0[geom0['Element'].str.contains("Li") == False]
    #geom1 = geom1[geom1['Element'].str.contains("Li") == False]
    #!! 

    #drop all columns except w,v,u
    geom0.drop(geom0.columns[[0,1]], axis=1, inplace=True)
    geom1.drop(geom1.columns[[0,1]], axis=1, inplace=True)
   
    #convert pandas to numpy arrays
    geom0, geom1  = geom0.values, geom1.values
   
    #absolute atom movement from geom0 to geom1
    geomdiff = np.absolute(geom0 - geom1)

    #-1 in reshape infers the lengths of the array
    #difference in w column only
    w_geomdiff = geomdiff[:,2].reshape(-1,2)
    
    #summing diff in H and diff in Li, total difference in w
    tot_w_geomdiff = np.sum(w_geomdiff,axis=1)
    
    #determine atoms in each layer
    atms_per_layer = (len(tot_w_geomdiff))/layerno
    
    #reshape array into each atomic layer in the w direction
    tot_w_geomdiff_lay = tot_w_geomdiff.reshape(-1,atms_per_layer)
   
    
    #start std_deviation
    for j in range(layerno):
        std_dev = np.std(tot_w_geomdiff_lay[j])
        #print std_dev as string as otherwise it is rounded
        std_dev_print = 'std_dev_layer%d = %s' %((j+1), std_dev)
        #print std_dev_print
    #end std_deviation
    

    #average movement in w per layer
    tot_w_geomdiff_avg = (np.sum(tot_w_geomdiff_lay,axis=1))/atms_per_layer
   
    
    #in Angstrom
    #c = ((layerno-1)*2.0042) + 12
    #print 'c= ' + str(c)
    #halfcell
    c = 22.0210306 
    #wholecell
    #c = 32.0420257
    #44Ang box
    #c = 44.0420257
    
    tot_w_geomdiff_avg = tot_w_geomdiff_avg * c
    
    #outputting to list for pyplot
    tot_w_geomdiff_avg = tot_w_geomdiff_avg.tolist() 
    
    return tot_w_geomdiff_avg

def element_atomic_disp_calc(layerno,Element):
   #! used to denoted how code is different to tot_atomic_disp_calc

    geom0 = pd.DataFrame((infile(layerno)),
            columns=['Element','ionNumber','u','v','w'])

    #change dtype of columns
    geom0[['ionNumber', 'u', 'v','w']] = geom0[['ionNumber','u','v','w']].astype(float)

    """ if the order of table is not H Li H Li etc.
    then you have to round to fewer d.p. """
    
    #must round w as starting positions not all the same
    geom0['w'] = geom0['w'].round(3)
    geom0['u'] = geom0['u'].round(4)
    geom0['v'] = geom0['v'].round(4)
    
    #order geom0 based on w then u then v
    geom0 = geom0.sort_values(['w','u','v'],
                              ascending=[True,True,True])

    #find index order for sorted geom0
    geom0index = geom0.index.tolist()
    
    #geom1 is compared to geom0
    geom1 = pd.DataFrame((outfile(layerno)),
            columns=['Element','ionNumber','u','v','w'])

    geom1[['ionNumber', 'u', 'v','w']] = geom1[['ionNumber','u','v','w']].astype(float)
    
    #redorder geom1 using geom0index
    geom1 = geom1.reindex(geom0index)
   
    #!!have to get geom0 again w/o rounding!!
    geom0 = pd.DataFrame((infile(layerno)),
            columns=['Element','ionNumber','u','v','w'])


    geom0[['ionNumber', 'u', 'v','w']] = geom0[['ionNumber','u','v','w']].astype(float)
    
    geom0 = geom0.reindex(geom0index)
    #end of !!getting geom0 again!!
    
    #!
    #select only the rows containing the desire element
    geom0 = geom0[geom0['Element'].str.contains(Element) == True]
    geom1 = geom1[geom1['Element'].str.contains(Element) == True]

    #drop all columns except w,v,u
    geom0.drop(geom0.columns[[0,1]], axis=1, inplace=True)
    geom1.drop(geom1.columns[[0,1]], axis=1, inplace=True)
   
    #convert pandas to numpy arrays
    geom0, geom1  = geom0.values, geom1.values
   
    #absolute atom movement from geom0 to geom1
    geomdiff = np.absolute(geom0 - geom1)
    
    
    #determine atoms in each layer
    atms_per_layer = (len(geomdiff))/layerno

    #!
    #difference in w column only
    w_geomdiff = geomdiff[:,2]


    #!
    #reshape array into each atomic layer in the w direction
    w_geomdiff_lay = w_geomdiff.reshape(-1,atms_per_layer)

    #!
    #average movement in w per layer for element
    w_geomdiff_avg = (np.sum(w_geomdiff_lay,axis=1))/atms_per_layer
   
    
    #convert fract coords to Angstrom
    c = ((layerno-1)*2.0042) + 12

    #!
    w_geomdiff_avg = w_geomdiff_avg * c
    
    #!
    #outputting to list for pyplot
    w_geomdiff_avg = w_geomdiff_avg.tolist() 
    
    return w_geomdiff_avg

def plot2d(inital_lay,final_lay,interval=1):
    import matplotlib.pyplot as plt
    
    colors = (
            '#ff0000','#ffbf00','#ffff00',
            '#80ff00','#000000','#00ffff',
            '#0000ff','#8000ff','#ff00ff',
            '#ff0000','#ffbf00','#ffff00',
            '#80ff00','#000000','#00ffff',
            '#0000ff','#8000ff','#ff00ff'
            )

    X = []

    for lay in range(inital_lay,final_lay+1,interval):
        #!
        #Y = atomic_disp_calc(lay,'Li')
        Y = LiH_rippling(lay)
        
        no_of_lay = (len(Y))
        if (no_of_lay % 2) == 0:
            #slice list to include half layers
            Y = Y[0:(no_of_lay/2)]

        else:
            #slice odd number of layer to include only unique items
            Y = Y[0:((no_of_lay+1)/2)]

        for i in range (len(Y)):
            X.append(i)
            
        plt.plot(X,Y,
                 linestyle='-',
                 marker='o',
                 ms=7,
                 color=colors[(lay-2)],
                 label = str(lay)+'L'
                )  
 
        #reset X for next plot
        X = []
    
    plt.grid(True)
    plt.xlabel('Layers from surface',fontsize=20)
    plt.ylabel('LiH Displacement /Ang',fontsize=20)
    plt.axis([-0.1,6.1,-0.001,0.08])#x_start,x_stop,y_start,y_stop
    plt.legend(loc=1,fontsize=20)
    plt.show()
    
    #custom_xticks = ['+5','+4','+3','+2','+1','0','-1','-2','-3','-4','-5']
    #plt.xticks(data[:,0],custom_xticks,fontsize=14)

def plot3darrows(layerno):
    from mayavi.mlab import *
    geom0 = pd.DataFrame((infile(layerno)),
                         columns=['Element','ionNumber','u','v','w'])
    geom1 = pd.DataFrame((outfile(layerno)),
                         columns=['Element','ionNumber','u','v','w'])
     
    geom0[['ionNumber', 'u', 'v','w']] = geom0[['ionNumber','u','v','w']].astype(float)
    
    geom1[['ionNumber', 'u', 'v','w']] = geom1[['ionNumber','u','v','w']].astype(float)
    
    element = ['H','Li','N']
    colors = [(1,1,1),(1,0,0),(0,0,1)]
    #seperate dataframe into H and Li
    for j in range(2):  
        #select single element
        geom0_el = geom0[geom0['Element'].str.contains(element[j]) == True]
        geom1_el = geom1[geom0['Element'].str.contains(element[j]) == True]

        #drop columns
        geom0_el.drop(geom0_el.columns[[0,1]], axis=1, inplace=True)
        geom1_el.drop(geom1_el.columns[[0,1]], axis=1, inplace=True)
        
        geomdiff_el = geom1_el - geom0_el #usually geom0_el - geom1_el
        #rename columns
        geomdiff_el.columns = ['du','dv','dw']

        quiver3d(geom0_el['u'],geom0_el['v'],geom0_el['w'],
                 geomdiff_el['du'],geomdiff_el['dv'],
                 geomdiff_el['dw'],mode='arrow',
                 color=colors[j],resolution=32)

    outline(color=(0,0,0))

    #adding fixed N atoms to 3d plot
    geom1_N = geom1[geom1['Element'].str.contains('N') == True]
    geom1_N = geom1_N.reset_index()
    for counter in range(2):
        points3d(geom1_N.loc[counter,'u'],
                 geom1_N.loc[counter,'v'],
                 geom1_N.loc[counter,'w'],
                 color=(0,0,1),resolution='32',
                 scale_factor=.02
                )
