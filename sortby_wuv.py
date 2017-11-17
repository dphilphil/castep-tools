import numpy as np
import pandas as pd

f = 'demo.xyz'

#grab data
raw = np.genfromtxt(f, delimiter=',')
unordered = pd.DataFrame((raw), columns=['Element','u','v','w'])
#order
unordered[['u', 'v','w']] = unordered[['u','v','w']].astype(float) #change dtype
ordered = unordered.sort_values(['w','u','v'], ascending=[True,True,True])
#export
ordered.to_csv(str(f)+".ordered", sep="\t",header=False)
