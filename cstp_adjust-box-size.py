import numpy as np

my_data = np.genfromtxt('LiHinput.cell', delimiter=',')
my_data[:,3] = (my_data[:,3]*32.0420257)/44.0420257
my_data[:,3] = my_data[:,3] + 0.1362335 #shift atoms so they're central in box
my_data[:,3] = np.round(my_data[:,3],7)
np.savetxt('LiHoutput.cell',my_data,delimiter='\t',fmt='%1.7f')

