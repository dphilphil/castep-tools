import numpy as np

#Ammonia Fractional Coordinates N followed by 3H
"""
OrigAmmonia = np.array([[0.5000000,0.5000000,0.8173434],
                        [0.5785366,0.5000000,0.8261927],
                        [0.4607241,0.5680054,0.8261958],
                        [0.4607241,0.5680054,0.8261958],
                       ])
"""
OrigAmmonia = np.array([[0.5000000,0.5000000,0.7948414],
                        [0.5785366,0.5000000,0.8036907],
                        [0.4607241,0.4319946,0.8036938],
                        [0.4607241,0.5680054,0.8036938],
                       ])

	
def Rotate_and_Invert(EndAngle,Step):
    #Setup blank RotatedAmmonia Array, getting Z axis as doesn't change
    RAmmonia = np.array([[0,0,OrigAmmonia[0,2]],
                         [0,0,OrigAmmonia[1,2]],
                         [0,0,OrigAmmonia[2,2]],
                         [0,0,OrigAmmonia[3,2]]
                        ])

    Counter = EndAngle/Step
    
    #make blank array with xyz
    Ammonia_All = np.zeros((4,3))

    for AngleMultiplier in range(Counter+1):
        Angle = AngleMultiplier*Step
        Angle = np.radians(Angle)
                
        for i in range(4):
            x, y  = OrigAmmonia[i,0], OrigAmmonia[i,1] #Z is fixed!
            
            #Origin pq
            p,q = 0.5, 0.5
            
            #benn.org/2007/01/06/
            Xprime = p + (np.cos(Angle)* (x-p)) - (np.sin(Angle) * (y-q))
            Yprime = q + (np.sin(Angle)* (x-p)) + (np.cos(Angle) * (y-q))
            
            Xprime,Yprime = np.round(Xprime,7), np.round(Yprime,7)
            
            RAmmonia[i,0] = Xprime
            RAmmonia[i,1] = Yprime
        Ammonia_All = np.append(Ammonia_All,RAmmonia,axis=0)
    
    Ammonia_All = Ammonia_All[4:]

    #inversion code
    Dist2Inversion = np.array(Ammonia_All - 0.5) #0.5 is inversion point as xyz =0.5,0.5,0.5
    InvertedAmmonia = 0.5 + (Dist2Inversion*-1)
    Checker = (InvertedAmmonia+Ammonia_All)/2 #should =0.5
    print "InternalCheck. All values should equal the Inversion Point: \n", Checker

    for AngleMultipler in range(Counter+1):
        mini = AngleMultipler*4
        maxi = (AngleMultipler+1)*4
        
        print str(AngleMultipler*Step) +' deg'
        print Ammonia_All[mini:maxi]
        print InvertedAmmonia[mini:maxi]
