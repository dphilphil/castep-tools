import numpy as np
import os

#usr input
castepf = 'GOPT_L-BFGS_LiH33_7LAYrepresenting13LAY_wFixed2BottomLays_750eV_kpn551_30AngBox_dipolecorr_NH3startingat2.3_000deg.castep'
li_atom = '59' #Li below ammonia
ab = 12.0251939
c = 30.0251951

#openfil and append data to lists
search = open(castepf, "r")

#empty lists
licoord = []
ncoord =[]
energies = []

for line in search:
    if "Final energy" in line:
        energy = line.split(' ') [17] # split on whitespaces and select energy from line
        energies.append(energy)
    
    if "x  Li          %s" % (li_atom)  in line:
        print (line.split(' '))
        licoord.append(line.split(' ') [33]) #append x
        licoord.append(line.split(' ') [36]) #   "   y
        licoord.append(line.split(' ') [39]) #   "   z

    if "x  N" in line:
        ncoord.append(line.split(' ') [35]) #append x 
        ncoord.append(line.split(' ') [38]) #   "   y 
        ncoord.append(line.split(' ') [41]) #   "   z

search.close()

#convert lists to arrays
energies = np.asarray(energies,dtype=float)

nitrogen = np.asarray(ncoord, dtype=float)
nitrogen = np.reshape(nitrogen,(-1,3))
#convert to absolute positions
nitrogen[:,:2]*= ab
nitrogen[:,2]*= c

lithium = np.asarray(licoord, dtype=float)
lithium = np.reshape(lithium,(-1,3))
lithium[:,:2]*= ab
lithium[:,2]*= c

#cal Li_H bondlength
lenbond  = np.sqrt(np.sum( ((nitrogen - lithium)**2), axis=1))

"""
#The last bond length is deleted as CASTEP outputs the coordinates of the last geometry step twice.
#Once as a geom step (with energy) and again as the 'Final Configuration' without energy.
"""
lenbond = np.delete(lenbond,-1,axis=0)

#outfil
outstack = np.column_stack((energies,lenbond))

basename = str(os.path.splitext(castepf)[0])
fname = 'Energy-v-lithium-ammonia-dist_' + basename + '.txt'

np.savetxt(fname,outstack,fmt='%s',header='Energy/eV \tLiN bondleng/Ang',delimiter='\t')
