import numpy as np
import matplotlib.pyplot as plt
import os

#IMPORTATION DES DONNES METEOS (VARIABLES EN FONCTION DU TEMPS)
fname = "meteo.csv"
f = open (fname)
data = f.read()
f.close()

lines = data.split('\n')
header = lines[0].split(';')
#lines = lines[1:]
print(header)

# float_data shape is (time,features)
float_data=np.zeros((len(lines),len(header)-1))
for i,line in enumerate(lines):
    values = [float(x) for x in line.split(';')[1:]]
    print(line.split(';')[0:])
    print(line)
    print(values)
    float_data[i,:]=values

#print (float_data[0,:])
print("data shape is {}".format(float_data.shape))
















              
       
        
    
        
        
        
        
        
     
        







