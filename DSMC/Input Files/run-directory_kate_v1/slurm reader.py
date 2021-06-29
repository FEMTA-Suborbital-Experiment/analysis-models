import pandas as pd
import matplotlib.pyplot as plt

#User Defined Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Slurm file version
version = '80558'

#Labels for stats_style NO SPACES
head=['Step', 'c_conv_sum', 'Np', 'Npave', 'Nexit', 'Nexitave', 'Nscoll', 'Nscollave', 'Nparent', 'Nchild', 'Nsplit']


#Timesteps per Output 
tpo = 1

#Convergence Criterion
error = 1E-4;

#Note: This code assumes grid refine happens every 100 timesteps
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''
DO NOT TOUCH ANYTHING BELOW THIS LINE
'''

header = ' '.join(head)
data = pd.DataFrame(columns = head)

with open('slurm-'+version+'.out') as f:
    content = f.readlines()
content = [x.strip() for x in content] 

dataStart = []
for i in range(len(content)):
    '''
    This for loop find the first starting point
    '''
    if content[i] == header:
        dataStart.append(i+1)
        break

for i in range(i+2,len(content)):
    if content[i] == header:
        dataStart.append(i+2)

for i in dataStart:
    if i == dataStart[0]:
        rownum=100/tpo+1  
    else:
        rownum=100/tpo 
    
    fileread = pd.read_csv('slurm-'+version+'.out',names=head,skiprows=i, nrows=rownum,delim_whitespace=True, header = None)
    data = data.append(fileread,ignore_index=True) 

ts = data['Step'].tolist()
#Np = data['Np'].tolist()
Np = data['c_conv_sum'].tolist()

#Remove Errors if applicable 
try: 
    int(Np[-1])
except ValueError:
    ts.pop(-1)
    Np.pop(-1)

ts = [int(i) for i in ts]
Np = [int(i) for i in Np]

#Plot Data
plt.plot(ts,Np)
plt.xlabel('Timestep')
plt.ylabel('Number of Particles')

#Test Convergence 
if abs(Np[-1]-Np[-2])/Np[-1]<=error:
    for i in range(1,len(ts)+1):
        if abs(Np[i]-Np[i-1])/Np[i]<=error:
            print('\nConverged at timestep '+str(ts[i]))
            
            plt.scatter(ts[i],Np[i],marker='*',color='r')
            plt.legend(['Np','Converged Timestep'])
            break
else:
    print('\nNot Converged :(\n')