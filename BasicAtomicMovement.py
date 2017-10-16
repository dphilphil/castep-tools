"""  
/   ******************************************************** 
    *                                                      *        
    * Prog:    BasicAtomicMovement		           *   
    * Author:  Phillip Marks                               *   
    * Date:    16-10-2017                                  *     
    * Source:  https://github.com/dphilphil/castepy/       *   
    *                                                      *   
    ********************************************************  /  
"""

import pandas as pd
import numpy as np

ab = 17.0743920664123
c =  35.0743920664123 

#import
inp = pd.read_csv('x.cell',\
sep=r"\s+", names=['atom','ind','x1','y1','z1',])

out = pd.read_csv('x-out.cell',\
sep=r"\s+", names=['atom','ind','x2','y2','z2'])

#multiply by abc
inp['x1'] = inp['x1'] * ab
inp['y1'] = inp['y1'] * ab
inp['z1'] = inp['z1'] * c
out['x2'] = out['x2'] * ab
out['y2'] = out['y2'] * ab
out['z2'] = out['z2'] * c

#declare x1 ... z2

x1 = inp['x1']
y1 = inp['y1']
z1 = inp['z1']
x2 = out['x2']
y2 = out['y2']
z2 = out['z2']

d = ((x2 - x1)**2) + ((y1-y2)**2) + ((z1-z2)**2)
d = np.sqrt(d)

print inp, out
print d
