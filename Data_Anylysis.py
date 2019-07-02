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
    
#sum_ALL_MI =[]
def time_shift(T):
    
    xsT=[]
    ysT=[]
    sum_ALL_MI =[]
    A32mBS20BH_ETH_T=[]
    mutu_info_ET_t=[]
    ysvT=[]
    ysv=[]
    for i in range (4025-T):
        A32mBS20BH_ETH_T.append(A32mBS20BH_ETH[i+T])
    
    
    
    for i in range (68):
        ysvT=[]
        xsT = np.asarray(A32mBS20BH_ETH_T) + 1
        ysv=np.swapaxes(A32mBS20BH_ACT[i:(i+1), :], 0,1)
        ysv=ysv.reshape(4026, )
        for i in range (4025-T):
            ysvT.append(ysv[i])
        ysT= np.asarray(ysvT) + 10
        
        
        k=pyinform.mutualinfo.mutual_info(xsT, ysT, bx=0, by=0, b=2.0, local=False)
        mutu_info_ET_t.append(k)
        
    we=sum(mutu_info_ET_t)
    return we
    #print(len(ysT),len (xsT))
    #plt.plot(sum_ALL_MI)
# Calling above Shift function
numset=[]
for i in range (100):
    numset.append(5*i)
sum_ALL_MI=[]
qaz=[]
for event in numset:
    qaz=time_shift(event)
    sum_ALL_MI.append(qaz)
    
plt.plot(sum_ALL_MI)


# Time Shift of T, Top 5 elements only
sum_ALL_MI =[]
def time_shift_5(T):
    
    xsT=[]
    ysT=[]
    TP_5_MI =[]
    A32mBS20BH_ETH_T=[]
    mutu_info_ET_t=[]
    ysvT=[]
    ysv=[]
    io=0
    for i in range (4026-T):
        A32mBS20BH_ETH_T.append(A32mBS20BH_ETH[i+T])
       

    ArT=[]
    ArryHT_E_Shift = np.asarray(ArT)
    for i in range (68):
        ysvT=[]
        xsT = np.asarray(A32mBS20BH_ETH_T) + 2
        ysv=np.swapaxes(A32mBS20BH_ACT[i:(i+1), :], 0,1)
        ysv=ysv.reshape(4026, )
        ut=i
        for i in range (4026-T):
            ysvT.append(ysv[i])
        ysT= np.asarray(ysvT) + 10
        u9=len(ysT)
        
        k=pyinform.mutualinfo.mutual_info(xsT, ysT, bx=0, by=0, b=2.0, local=False)
        mutu_info_ET_t.append(k)
        
    mutu_info_ET_t.sort(reverse = True )
    TP_5_MI=mutu_info_ET_t[:5]
    we=np.asarray(TP_5_MI)
    we=we.reshape(1,-1)
    
    
    #return ArryHT_E_Shift.shape, we.shape
    return we


#Heat map plot

plt.rcParams["figure.figsize"] = 10,10
fig, ax = plt.subplots()
PLT_HT=pd.DataFrame(AWRY)
sns.heatmap(PLT_HT, cmap="plasma")


##Saving the Heatmap

o=sns.heatmap(PLT_HT, cmap="plasma")
figure = o.get_figure()    
figure.savefig('ETH_1000_shift(top_5)_heatmap.png', dpi=400)
