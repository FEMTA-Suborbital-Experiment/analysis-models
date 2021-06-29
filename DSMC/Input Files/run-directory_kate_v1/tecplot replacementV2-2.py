'''
NOTE: this code is for axisymetric only

NOTE FEMTA POST PROCESSING IS BROKEN WILL FIX LATER - Cathode is fine
'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time
from numba import jit


start_time = time.time()

#User Defined Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ts = '2400' # Timestep

'''
Make sure you download the corresponding file
    1 = Yes
    0 = No
'''
femtaSwitch = 0 # femtaData.out  <<<<< Broken for now do not turn to 1
cathodeSwitch = 1 # cathodeData.out
plotSwitch = 1 # Would you like a semi-nice looking contour plot?

# Cool Numbers
cathodeDist = 0.005 # Distance from FEMTA surface (no Cavity) to Cathode Beginning (m)

# Fixed Numbers
radius = 0.0127 # Radius of Cathode (m)
cavityDepth = 0.00324 # Cavity Depth of FEMTA (m)
cathodeLength = 0.0281 # Length of Cathode (m)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''
DO NOT TOUCH ANYTHING BELOW THIS LINE
'''

#Cathode Geometry & General Info
cathigh = -(cathodeDist + cavityDepth)
catlow = -(cathodeDist + cavityDepth + cathodeLength)

femtaHead = ['xc', 'yc', 'xlo', 'ylo', 'zlo', 'xhi', 'yhi', 'zhi', 'f_femtaPress_fix[1]', 'f_femtaPress_fix[2]', 'f_femtaNrho_fix[1]', 'f_femtaNrho_fix[2]', 'f_femtaNrho_fix[3]', 'f_femtaNrho_fix[4]', 'f_femtaNrho_fix[5]']
cathodeHead = ['xc', 'yc', 'f_cathodePress_fix[1]', 'f_cathodePress_fix[2]', 'f_cathodeNrho_fix[1]', 'f_cathodeNrho_fix[2]', 'f_cathodeNrho_fix[3]', 'f_cathodeNrho_fix[4]']

'''
Press_fix:
    fix[1]: Pressure
    fix[2]: Temperature

Nrho_fix
    fix[1]: Particle Count
    fix[2]: x-velocity
    fix[3]: y-velocity
    fix[4]: Number Density
'''

#Numba Jit (Pandas Excluded because not compatible with jit)
@jit(nopython=True)
def computeCathode(xc,yc,press,temp,dataSize):
# Computes pressure and temperature inside the physical cathode only!    
        
    cathodePress = []
    cathodeTemp = []
    
    cathX = []
    cathY = []
    
    for i in range(dataSize):
        x = xc[i]
        y = yc[i]
    
        if x >= catlow and x <= cathigh and abs(y) <= radius:
            cathodePress.append(press[i])
            cathodeTemp.append(temp[i])
            cathX.append(x)
            cathY.append(y)

#        if i%10000 == 0:      
            #print('Cathode: {:.2f}% Complete'.format(i/dataSize*100))
    
    return(cathodePress,cathodeTemp,cathX,cathY)

#Switch Statements
if femtaSwitch == 1:
    femtaData = pd.DataFrame(columns = femtaHead)
    fileread = pd.read_csv('femtaData.'+ts+'.out',skiprows=9,names=femtaHead,delim_whitespace=True, header = None)
    femtaData = femtaData.append(fileread,ignore_index=True)
    
    femtaPress = femtaData['f_femtaPress_fix[1]'].to_numpy()
    femtaTemp = femtaData['f_femtaPress_fix[2]'].to_numpy()

if cathodeSwitch == 1:
    cathodeData = pd.DataFrame(columns = cathodeHead)
    fileread = pd.read_csv('cathodeData.'+ts+'.out',skiprows=9,names=cathodeHead,delim_whitespace=True, header = None)
    cathodeData = cathodeData.append(fileread,ignore_index=True)

    #Extracting useful data to lists
    dataSize = len(cathodeData)
    xc = cathodeData['xc'].to_numpy()
    yc = cathodeData['yc'].to_numpy()
    cPress = cathodeData['f_cathodePress_fix[1]'].to_numpy()
    cTemp = cathodeData['f_cathodePress_fix[2]'].to_numpy()

    (cathodePress,cathodeTemp,cathX,cathY) = computeCathode(xc,yc,cPress,cTemp,dataSize)
    
if plotSwitch == 1:
    if cathodeSwitch != 1:
        cathodeData = pd.DataFrame(columns = cathodeHead)
        fileread = pd.read_csv('cathodeData.'+ts+'.out',skiprows=9,names=cathodeHead,delim_whitespace=True, header = None)
        cathodeData = cathodeData.append(fileread,ignore_index=True)
        
        dataSize = len(cathodeData)
        xc = cathodeData['xc'].to_numpy()
        yc = cathodeData['yc'].to_numpy()
        cPress = cathodeData['f_cathodePress_fix[1]'].to_numpy()
        cTemp = cathodeData['f_cathodePress_fix[2]'].to_numpy()
        
        
    (cathodePress,cathodeTemp,cathX,cathY) = computeCathode(xc,yc,cPress,cTemp,dataSize)
    
    fig=plt.figure(1)
    plt.tricontourf(xc,yc, cPress,50)
    plt.title('Pressure (Pa)') 
    plt.axvline(x=catlow,c='r')
    plt.axvline(x=cathigh,c='r')
    plt.colorbar()
    
    plt.figure(2)
    plt.tricontourf(xc,yc, cTemp,50)
    plt.title('Temperature (K)')
    plt.axvline(x=catlow,c='r')
    plt.axvline(x=cathigh,c='r')
    plt.colorbar()

#OUTPUTS
if femtaSwitch == 1:
    print('\n----FEMTA RESULTS----\n')
    print('Max Femta Pressure: '+str(np.percentile((femtaPress,95))))
    print('Max Femta Temperature: '+str(np.percentile((femtaTemp,95))))
if cathodeSwitch == 1:
    print('\n----CATHODE RESULTS----\n')   
    print('95th Percentile Cathode Pressure: '+str(np.percentile(cathodePress,95)))
    print('Average Cathode Pressure: '+str(np.mean(cathodePress)))
    print('Max Cathode Temperature: '+str(np.percentile(cathodeTemp,95)))

print("\nTime Elapsed: %.4f seconds" % (time.time() - start_time))