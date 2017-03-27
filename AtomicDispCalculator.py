import pandas as pd
import numpy as np
import linecache


def infile(layerno):
   
    castepf = 'GOPT_L-BFGS_LiH33_11LAY_750eV_kpn551_wFixedBulk_44AngBox_SymmGen_w2NH3at_3AngfromSurf_150deg.castep'

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

def outfile(layerno):

    castepf = 'GOPT_L-BFGS_LiH33_11LAY_750eV_kpn551_wFixedBulk_44AngBox_SymmGen_w2NH3at_3AngfromSurf_150deg.castep'

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

#calculate total atomic displacement in layers
def tot_atomic_disp_calc(layerno,Ammonia=True):
    if Ammonia==True: print "add number of NH3 atoms(6) to layers?"
    
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
    c = 44.0420257
    
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
        Y = tot_atomic_disp_calc(lay)
        
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
    
def alternativeplot2d(inital_lay,final_lay,interval=1):
    
    colors = (
            '#ff0000','#ffbf00','#ffff00',
            '#80ff00','#000000','#00ffff',
            '#0000ff','#8000ff','#ff00ff',
            '#ff0000','#ffbf00','#ffff00',
            '#80ff00','#000000','#00ffff',
            '#0000ff','#8000ff','#ff00ff'
            )
             
    for lay in range(inital_lay,final_lay +1,interval):
        Y = tot_atomic_disp_calc(lay)

        no_of_lay = len(Y)
        
        #calculate X
        Xstart = (-no_of_lay/2.0) + 0.5
        Xend = (no_of_lay/2.0) -0.5
        #numpy array used as X axis has floats so can't use range()
        X = np.arange(Xstart,Xend+1,1)
        
        plt.plot(X,np.log10(Y),
                linestyle='-',
                marker='o',
                ms=7,
                color=colors[(lay-2)],
                label = str(lay) +'L'
                )
        
    plt.grid(True)
    plt.xlabel('Number of layers from centre',fontsize=16)
    plt.ylabel('log10(LiHDeltas)',fontsize=16)
    plt.axis([-5,5,-5.5,-1.0])
    plt.legend(loc=4, fontsize=16)
    plt.show()

    #custom_xticks = ['+5','+4','+3','+2','+1','0','-1','-2','-3','-4','-5']
    #plt.xticks(data[:,0],custom_xticks,fontsize=14)
