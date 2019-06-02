import pandas as pd 
import h5py
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pyinform.dist import Dist
from pyinform.shannon import mutual_info
import random
import pyinform.mutualinfo 


A32m=h5py.File('AML32_moving.hdf5', 'r')
A32m.keys()
A32mBS20= A32m['BrainScanner20170424_105620']
A32mBS20BH= A32mBS20['Behavior']
A32mBS20BH_ETH= np.array(A32mBS20BH.get('Ethogram'))
A32mBS20BH_CMV= np.array(A32mBS20BH.get('CMSVelocity'))
A32mBS20BH_AV= np.array(A32mBS20BH.get('AngleVelocity'))

plt.plot(A32mBS20BH_ETH)
plt.rcParams["figure.figsize"] = 10,10
fig, ax = plt.subplots()
sns.heatmap([A32mBS20BH_ETH])

A32mBS20NR= A32mBS20['Neurons']
A32mBS20BH_ACT= np.array(A32mBS20NR.get('Activity'))

#Heat map plot
A32mBS20BH_ACT_plt=pd.DataFrame(A32mBS20BH_ACT)
sns.heatmap(A32mBS20BH_ACT_plt, cmap="Greens")

# Select multiple rows from index x to z
A32mBS20BH_ACT_rows_30_35 = A32mBS20BH_ACT[30:31, :]

A32mBS20BH_ACT_rows_30_35_plt=pd.DataFrame(A32mBS20BH_ACT_rows_30_35)
sns.heatmap(A32mBS20BH_ACT_rows_30_35, cmap="Greens")

#Plot two histograms at the same time 

xs=A32mBS20BH_ETH
# for i in range (68):
    
#     ys=np.swapaxes(A32mBS20BH_ACT[i:(i+1), :], 0,1)
#     ys=ys.reshape(4026, )

#     bins = np.linspace(-10, 10, 20)

#     plt.figure()
#     plt.hist(xs, bins, alpha=0.5, label='x')
#     plt.hist(ys, bins, alpha=0.5, label='y')
#     plt.legend(loc='upper right')
#         #pyplot.show()
#     plt.savefig("graph"+str(i) + '.png')
#     plt.close()


# MUTUAL info
#pyinform.mutualinfo.mutual_info(xs, ys, bx=0, by=0, b=2.0, local=False)


# Mutual Info for Activity array and Ethogram
mutu_info=[]
for i in range (68):
    ys=np.swapaxes(A32mBS20BH_ACT[i:(i+1), :], 0,1)
    ys=ys.reshape(4026, )
    ysP= ys+ 10
    k=pyinform.mutualinfo.mutual_info(xsP, ysP, bx=0, by=0, b=2.0, local=False)
    mutu_info.append(k)
    
#Mutual Info for Activity array and CMS Velocity     
xsV=A32mBS20BH_CMV + 1
mutu_info_CMSV=[]
for i in range (68):
    ysv=np.swapaxes(A32mBS20BH_ACT[i:(i+1), :], 0,1)
    ysv=ysv.reshape(4026, )
    ysV= ysv+ 10
    k=pyinform.mutualinfo.mutual_info(xsV, ysV, bx=0, by=0, b=2.0, local=False)
    mutu_info_CMSV.append(k)

#Mutual Info for Activity array and Angular Velocity 
xsAV=A32mBS20BH_AV + 3
mutu_info_AV=[]
for i in range (68):
    ysv=np.swapaxes(A32mBS20BH_ACT[i:(i+1), :], 0,1)
    ysv=ysv.reshape(4026, )
    ysAV= ysv+ 10
    k=pyinform.mutualinfo.mutual_info(xsAV, ysAV, bx=0, by=0, b=2.0, local=False)
    mutu_info_AV.append(k)
    
plt.plot(mutu_info_AV)


## making the time shift of T in Etho. data
A32mBS20BH_ETH_T50=[]

for i in range (3975):
    A32mBS20BH_ETH_T50.append(A32mBS20BH_ETH[i+50])
    
#Mutual Info with time shift

xsT50 = np.asarray(A32mBS20BH_ETH_T50) + 1

mutu_info_ET_t50=[]
ysvT50=[]
ysv=[]
for i in range (68):
    ysvT50=[]
    ysv=np.swapaxes(A32mBS20BH_ACT[i:(i+1), :], 0,1)
    ysv=ysv.reshape(4026, )
    for i in range (3975):
        ysvT50.append(ysv[i])
    
    ysT50= np.asarray(ysvT50) + 10
    k=pyinform.mutualinfo.mutual_info(xsT50, ysT50, bx=0, by=0, b=2.0, local=False)
    mutu_info_ET_t50.append(k)

