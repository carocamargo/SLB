#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 19:22:23 2022

@author: ccamargo
"""

import xarray as xr
import numpy as np

#%% get budget components

period = ['1993-2017'] # full years
y0,y1=period[0].split('-')
t0='{}-01-01'.format(int(y0))
t1='{}-12-31'.format(int(y1)-1)
path = '/Volumes/LaCie_NIOZ/data/budget/trends/' 
    
#%% alt
comp = 'alt'
file = comp+'.nc'
ds=xr.open_dataset(path+file)
ds = ds.sel(periods=period)
idx = np.where(ds.ICs =='bic_tp')[0][0]
alt_trends = np.array(ds.best_trend[:,idx,:,:])
alt_uncs = np.array(ds.best_unc[:,idx,:,:])
alt_names = ['ENS']

#%% dyn
comp = 'dynamic'
file = comp+'.nc'
ds=xr.open_dataset(path+file)
ds = ds.sel(periods=period)
idx = np.where(ds.ICs =='bic_tp')[0][0]
dyn_trends = np.array(ds.best_trend[:,idx,:,:])
dyn_uncs = np.array(ds.best_unc[:,idx,:,:])
dyn_names=['ENS']
#%% steric
comp = 'steric'
file = comp+'.nc'
ds=xr.open_dataset(path+file)
ds = ds.sel(periods=period)
ste_trends = np.array(ds.trend_full[:,0,:,:])
ste_uncs = np.array(ds.unc[:,0,:,:])
ste_names = np.array(ds.names)
#%% barystatic
comp = 'barystatic'
file = comp+'.nc'
ds=xr.open_dataset(path+file)
ds = ds.sel(periods=period)
ds = ds.where((ds.lat>-66) & (ds.lat<66),np.nan)
bar_trends = np.array(ds.SLA[:,0,:,:])
bar_uncs = np.array(ds.SLA_UNC[:,0,:,:])      
bar_names = np.array(ds.names)

#%%
n_pos = alt_trends.shape[0] * ste_trends.shape[0] * bar_trends.shape[0] * dyn_trends.shape[0]
res = np.zeros((n_pos,alt_trends.shape[1],alt_trends.shape[2]))
unc = np.full_like(res,0)

ipos=0
combs=[]
names = np.full_like(np.zeros((n_pos,4)),np.nan).astype(str)
for ialt in range(alt_trends.shape[0]):
    for iste in range(ste_trends.shape[0]):
        for ibar in range(bar_trends.shape[0]):
            for idyn in range(dyn_trends.shape[0]):
                combs.append('{}_{}_{}_{}'.format(ialt,iste,ibar,idyn))
                names[ipos,:] = [alt_names[ialt], ste_names[iste], bar_names[ibar], dyn_names[idyn]]
                
                res[ipos] = np.array(alt_trends[ialt] - 
                                (ste_trends[iste] + 
                                bar_trends[ibar] + 
                                dyn_trends[idyn] )
                                                    )
                unc[ipos] = np.array(
                    np.sqrt(
                    np.abs(
                                alt_uncs[ialt]**2 - 
                                (ste_uncs[iste]**2 + 
                                bar_uncs[ibar]**2 + 
                                dyn_uncs[idyn]**2 ))
                                                    ))
                ipos=ipos+1

da = xr.Dataset(data_vars = {'res':(('comb','lat','lon'),res),
                             'unc':(('comb','lat','lon'),unc),
                             'names':(('combs','dataset'),names)
                             },
                coords={'comb':combs,
                        'dataset':['alt','ste','bar','dyn'],
                             'lat':ds.lat,
                             'lon':ds.lon}
                )

da.to_netcdf('/Volumes/LaCie_NIOZ/data/budget/combinations.nc')
