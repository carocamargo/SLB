#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 15:22:02 2022

@author: ccamargo
"""


import numpy as np
# import xarray as xr
# # import pickle
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.append("/Users/ccamargo/Documents/github/SLB/")

from utils_SLB import unc_test, agree_test, zeta_test
#%%
path = '/Volumes/LaCie_NIOZ/data/budget/'
dic = pd.read_pickle(path+'budget.pkl')
#%%
key ='som' # cluster
res = np.array(dic[key]['res']['trend'])# cluster residual trend
unc = np.array(dic[key]['res']['unc']) # cluster uncertainty
mask = np.array(dic[key]['mask']) # clusters mask
n = dic[key]['n']
df = dic[key]['df']
alt = np.array(dic['alt']['ts'])
comp = np.array(dic['sum']['ts'])
res = np.array(alt-comp)

#%%
fig = plt.figure(figsize=(20,10))
for i in range(n): 
    icluster = int(i+1)
    mask_tmp = np.array(mask)
    mask_tmp[np.where(mask_tmp!=icluster)] = np.nan
    mask_tmp[np.isfinite(mask_tmp)]= 1
    
    y_alt = np.nanmean(alt*mask_tmp,axis=(1,2)) - np.nanmean(alt*mask_tmp)
    y_comp = np.nanmean(comp*mask_tmp,axis=(1,2)) - np.nanmean(comp*mask_tmp)
    
    ax = plt.subplot(6,3,icluster)
    plt.plot(y_alt,label = 'alt',linewidth=2)
    plt.plot(y_comp, label='sum',linewidth=2,linestyle='--')
    plt.plot(y_alt-y_comp, label='res',alpha = 0.5)
    plt.title('Cluster {}'.format(icluster))
    plt.legend()
plt.tight_layout()
plt.show()
#%%
plt.pcolor(dic['alt']['trend']-dic['sum']['trend'],
           vmin=-1,vmax=1,cmap='RdBu')
plt.colorbar()
plt.show()
#%%
plt.pcolor(dic['som']['res']['trend'],
           vmin=-1,vmax=1,cmap='RdBu')
plt.colorbar()
plt.show()
#%%
plt.figure()
y = np.array(dic['alt']['trend']).flatten()
plt.hist(y,**kwargs,label='alt')
y = np.array(dic['sum']['trend']).flatten()
plt.hist(y,**kwargs,label='sum')
plt.legend()
plt.show()
#%%
plt.figure()
y = np.array(dic['som']['alt']['trend']).flatten()
plt.hist(y,**kwargs,label='alt')
y = np.array(dic['som']['sum']['trend']).flatten()
plt.hist(y,**kwargs,label='sum')
#%%
y = np.array(dic['alt']['trend']).flatten()
x = np.array(dic['sum']['trend']).flatten()
plt.scatter(x,y)
plt.xlabel('altimetry')
plt.ylabel('sum')
plt.show()
#%%
from scipy.stats import gaussian_kde
#%% scatter ts 
plt.figure(figsize=(15,15))
nrow=2
ncol=2

mask_tmp = np.array(mask)
mask_tmp[np.isfinite(mask_tmp)]=1

plt.subplot(nrow,ncol,1)
plt.title('1 degree')
y = np.array(dic['alt']['ts'] * mask_tmp).flatten()
x = np.array(dic['sum']['ts'] * mask_tmp).flatten() 
x  = np.array(x - np.nanmin(x))/ (np.nanmax(x) - np.nanmin(x))
y  = np.array(y - np.nanmin(y))/ (np.nanmax(y) - np.nanmin(y))

x = x[np.isfinite(y)]
y = y[np.isfinite(y)]
y = y[np.isfinite(x)]
x = x[np.isfinite(x)]
plt.scatter(x,y,alpha=0.5)
plt.plot([0, 1], [0, 1], ls="--", c="black",alpha=0.8)
plt.xlim([-0.1,1.1])
plt.ylim([-0.1,1.1])

##
plt.subplot(nrow,ncol,2)
plt.title('GMSL')
y = np.nanmean(np.array(dic['alt']['ts']) *mask_tmp,axis=(1,2)) 
x = np.nanmean(np.array(dic['sum']['ts']) *mask_tmp,axis=(1,2)) 
x  = np.array(x - np.nanmin(x))/ (np.nanmax(x) - np.nanmin(x))
y  = np.array(y - np.nanmin(y))/ (np.nanmax(y) - np.nanmin(y))

x = x[np.isfinite(y)]
y = y[np.isfinite(y)]
y = y[np.isfinite(x)]
x = x[np.isfinite(x)]
plt.scatter(x,y,alpha=0.5)

plt.plot([0, 1], [0, 1], ls="--", c="black",alpha=0.8)
plt.xlim([-0.1,1.1])
plt.ylim([-0.1,1.1])

##
plt.subplot(nrow,ncol,3)
plt.title('SOM')
for i in range(n): 
    icluster = int(i+1)
    mask_tmp = np.array(mask)
    mask_tmp[np.where(mask_tmp!=icluster)] = np.nan
    mask_tmp[np.isfinite(mask_tmp)]= 1
    
    y = np.nanmean(alt*mask_tmp,axis=(1,2)) 
    y  = np.array(y - np.nanmin(y))/ (np.nanmax(y) - np.nanmin(y))
    x = np.nanmean(comp*mask_tmp,axis=(1,2)) #- np.nanmean(comp*mask_tmp)
    x  = np.array(x - np.nanmin(x))/ (np.nanmax(x) - np.nanmin(x))

    x = x[np.isfinite(y)]
    y = y[np.isfinite(y)]
    y = y[np.isfinite(x)]
    x = x[np.isfinite(x)]

    plt.scatter(x,y,alpha=0.5)

plt.plot([0, 1], [0, 1], ls="--", c="black",alpha=0.8)
plt.xlim([-0.1,1.1])
plt.ylim([-0.1,1.1])

##
plt.subplot(nrow,ncol,4)
plt.title('delta Maps')
for i in range(dic['dmap']['n']): 
    icluster = int(i+1)
    mask_tmp = np.array(dic['dmap']['mask'])
    mask_tmp[np.where(mask_tmp!=icluster)] = np.nan
    mask_tmp[np.isfinite(mask_tmp)]= 1
    
    y = np.nanmean(alt*mask_tmp,axis=(1,2)) 
    y  = np.array(y - np.nanmin(y))/ (np.nanmax(y) - np.nanmin(y))
    x = np.nanmean(comp*mask_tmp,axis=(1,2)) #- np.nanmean(comp*mask_tmp)
    x  = np.array(x - np.nanmin(x))/ (np.nanmax(x) - np.nanmin(x))

    x = x[np.isfinite(y)]
    y = y[np.isfinite(y)]
    y = y[np.isfinite(x)]
    x = x[np.isfinite(x)]

    plt.scatter(x,y,alpha=0.5)
plt.plot([0, 1], [0, 1], ls="--", c="black",alpha=0.8)
plt.xlim([-0.1,1.1])
plt.ylim([-0.1,1.1])

plt.tight_layout()
plt.show()
#%%

#%% scatter trend
plt.figure(figsize=(15,5))
nrow=1
ncol=3

mask_tmp = np.array(mask)
mask_tmp[np.isfinite(mask_tmp)]=1

plt.subplot(nrow,ncol,1)
plt.title('1 degree')
y = np.array(dic['alt']['trend'] * mask_tmp).flatten() 
x = np.array(dic['sum']['trend'] * mask_tmp).flatten()
x  = np.array(x - np.nanmin(x))/ (np.nanmax(x) - np.nanmin(x))
y  = np.array(y - np.nanmin(y))/ (np.nanmax(y) - np.nanmin(y))

x = x[np.isfinite(y)]
y = y[np.isfinite(y)]
y = y[np.isfinite(x)]
x = x[np.isfinite(x)]

plt.scatter(x,y)
plt.xlim([-0.1,1.1])
plt.ylim([-0.1,1.1])
plt.plot([0, 1], [0, 1], ls="--", c=".3")

plt.subplot(nrow,ncol,2)
plt.title('SOM')
y = np.array(dic['som']['alt']['trend']).flatten()
x = np.array(dic['som']['sum']['trend']).flatten()
x  = np.array(x - np.nanmin(x))/ (np.nanmax(x) - np.nanmin(x))
y  = np.array(y - np.nanmin(y))/ (np.nanmax(y) - np.nanmin(y))

x = x[np.isfinite(y)]
y = y[np.isfinite(y)]
y = y[np.isfinite(x)]
x = x[np.isfinite(x)]

plt.scatter(x,y)
plt.xlim([-0.1,1.1])
plt.ylim([-0.1,1.1])
plt.plot([0, 1], [0, 1], ls="--", c=".3")

ax = plt.subplot(nrow,ncol,3)
plt.title('dMAP')
y = np.array(dic['dmap']['alt']['trend']).flatten()
x = np.array(dic['dmap']['sum']['trend']).flatten()

# y = np.array(dic['dmap']['df']['alt_tr'])
# x = np.array(dic['dmap']['df']['sum_tr'])

x = x[np.isfinite(y)]
y = y[np.isfinite(y)]

y = y[np.isfinite(x)]
x = x[np.isfinite(x)]

x  = np.array(x - np.nanmin(x))/ (np.nanmax(x) - np.nanmin(x))
y  = np.array(y - np.nanmin(y))/ (np.nanmax(y) - np.nanmin(y))
plt.scatter(x,y)
plt.plot([0, 1], [0, 1], ls="--", c=".3")
plt.xlim([-0.1,1.1])
plt.ylim([-0.1,1.1])

plt.tight_layout()
plt.show()


#%%
plt.subplot(nrow,ncol,2)
for i in range(n): 
    icluster = int(i+1)
    mask_tmp = np.array(mask)
    mask_tmp[np.where(mask_tmp!=icluster)] = np.nan
    mask_tmp[np.isfinite(mask_tmp)]= 1
    
    y = np.nanmean(alt*mask_tmp,axis=(1,2)) - np.nanmean(alt*mask_tmp)
    x = np.nanmean(comp*mask_tmp,axis=(1,2)) - np.nanmean(comp*mask_tmp)
    # y_res = np.array(y_alt - y_comp)
    plt.scatter(x,y)

plt.title('SOM')
plt.tight_layout()
ptlt.show()
#%%
y = np.array(dic['res']['trend']).flatten()
plt.hist(y,**kwargs,label='res 1 deg')

y = np.array(df['res_tr']).flatten()
plt.hist(y,**kwargs,label='res som')

plt.xlim((-5,5))
plt.legend()
plt.show()
#%% histogram gmsl
plt.figure()
mask_tmp = np.array(mask)
mask_tmp[np.isfinite(mask_tmp)]=1
y_alt = np.nanmean(alt*mask_tmp,axis=(1,2)) - np.nanmean(alt*mask_tmp)
y_comp = np.nanmean(comp*mask_tmp,axis=(1,2)) - np.nanmean(comp*mask_tmp)
y_res = np.array(y_alt - y_comp)

plt.subplot(2,1,1)
_,_,_ = plt.hist(y_alt,**kwargs,label='alt')
_,_,_ = plt.hist(y_comp,**kwargs,label='sum')
plt.legend()

plt.subplot(2,1,2)
plt.plot(y_alt,label = 'alt',linewidth=2)
plt.plot(y_comp, label='sum',linewidth=2,linestyle='--')
plt.plot(y_alt-y_comp, label='res',alpha = 0.5)
plt.show()
#%% hitogram clusters
fig = plt.figure(figsize=(20,10))
for i in range(n): 
    icluster = int(i+1)
    mask_tmp = np.array(mask)
    mask_tmp[np.where(mask_tmp!=icluster)] = np.nan
    mask_tmp[np.isfinite(mask_tmp)]= 1
    
    y_alt = np.nanmean(alt*mask_tmp,axis=(1,2)) - np.nanmean(alt*mask_tmp)
    y_comp = np.nanmean(comp*mask_tmp,axis=(1,2)) - np.nanmean(comp*mask_tmp)
    y_res = np.array(y_alt - y_comp)
    ax = plt.subplot(3,6,icluster)
    # plt.plot(y_alt,label = 'alt',linewidth=2)
    # plt.plot(y_comp, label='sum',linewidth=2,linestyle='--')
    # plt.plot(y_alt-y_comp, label='res',alpha = 0.5)
    plt.title('Cluster {}'.format(icluster))
    

    kwargs = dict(histtype='stepfilled', 
                  alpha=0.3, density=True, 
                  bins=40, ec="k")
    _,_,_ = plt.hist(y_alt,**kwargs,label='alt')
    _,_,_ = plt.hist(y_comp,**kwargs,label='sum')
    # _,_,_ = plt.hist(y_res,**kwargs)
    if icluster==18:
        plt.legend()
    


    # plt.legend()
plt.tight_layout()
plt.show()
#%%
    
    
n_cluster=np.unique(mask_comb[np.isfinite(mask_comb)])
    modos = np.zeros((len(n_cluster),len(da.time)))
    ts=np.zeros((len(var),len(da.time)))
    tr=np.zeros((len(var)))
    plt.figure(figsize=(20,10))
    budget_cluster=np.full_like(mask,np.nan)
    budget_res = np.full_like(mask,np.nan)
    budget_noise = np.full_like(mask,np.nan)
    
    for i in range(0,int(np.nanmax(mask_comb))):
        icluster = i+1
        mask=np.array(mask_comb)
        mask[np.where(mask!=icluster)]=np.nan
        mask[np.isfinite(mask)]=1
        # plt.pcolor(mask)
        # plt.show()
        ampl = np.abs(np.nanmax(da.sla[:,:,:,0]*mask) - np.nanmin(da.sla[:,:,:,0]*mask))
        plt.subplot(y,x,icluster)
        # plt.subplot(111)
        
        for ivar in range(len(var)):
            ts[ivar,:] = cluster_mean(np.array(da.sla[:,:,:,ivar]), mask,time=da.time,lat=da.lat,lon=da.lon)
            out = sl.get_OLS_trend(tdec,ts[ivar])
            tr[ivar]=out[0]
            plt.plot(ts[ivar,:],label='{}: {} '.format(var[ivar],np.round(tr[ivar],3)))
        plt.legend()
    #% %
        budget_cluster[np.where(mask==1)] = np.array(tr[0] - (tr[1]+tr[2]))
        budget_res[np.where(mask==1)] = (np.array(tr[0] - (tr[1]+tr[2])))/ampl
        budget_noise[np.where(mask==1)]= np.array(np.array(tr[0] - (tr[1]+tr[2]))/tr[0])
        plt.title('Cluster {}'.format(icluster))
    plt.show()