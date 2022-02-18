#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 09:49:13 2022

@author: ccamargo
"""

import xarray as xr
import numpy as np
import os
import pandas as pd
import sys
sys.path.append("/Users/ccamargo/Documents/github/SLB/")

from utils_SLB import cluster_mean, plot_map_subplots, sum_linear, sum_square, get_dectime
from utils_SLB import plot_map2 as plot_map

sys.path.append("/Users/ccamargo/Documents/py_scripts/")
import utils_SL as sl

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean as cm
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
cmap_trend = cm.cm.balance
cmap_unc = cm.tools.crop(cmap_trend,0,3,0)

#%% get budget components

dic = {}

period = ['1993-2017'] # full years
y0,y1=period[0].split('-')
t0='{}-01-01'.format(int(y0))
t1='{}-12-31'.format(int(y1)-1)
path = '/Volumes/LaCie_NIOZ/data/budget/trends/' 
path2 = '/Volumes/LaCie_NIOZ/data/budget/ts/' 
flist = [file for file in os.listdir(path) if not file.startswith('.')]

for file in flist:
    comp = file.split('.')[0]
    # print(comp)
    
    ds=xr.open_dataset(path+file)
    # print(ds)
    da = xr.open_dataset(path2+file)
    # print(da)
    ds = ds.sel(periods=period)
    da = da.sel(time=slice(t0,t1))

    if comp =='dynamic' or comp=='alt':
        idx = np.where(ds.ICs =='bic_tp')[0][0]
        trend = np.array(ds.best_trend[0,idx,:,:])
        unc = np.array(ds.best_unc[0,idx,:,:])
        
        if comp=='alt':
            ts = np.array(da.sel(names='ENS')['SLA']) # mm
            tdec = get_dectime(da.time)
            #% % land mask
            da = da['sla_ens'][0,:,:]
            da = da.where((ds.lat>-66) & (ds.lat<66),np.nan)
            # da.plot()
            landmask = np.array(da.data)
            landmask[np.isfinite(landmask)]=1

        else:
            ts = np.array(da['ens_dyn_v1'] * 1000) # mm
        

    elif comp=='steric':
        
        ds=ds.sel(names=['ENS'])
        trend = np.array(ds.trend_full[0,0,:,:])
        unc = np.array(ds.unc[0,0,:,:])
        ts = np.array(da.sel(names='ENS')['steric_full'] * 1000) # mm
        
        trend_up = np.array(ds.trend_up[0,0,:,:])
        unc_up = np.array(ds.unc[0,0,:,:])
        ts_up = np.array(da.sel(names='ENS')['steric_up']*1000)
        dic['steric_up'] = {'trend':trend_up,
                            'unc':unc_up,
                            'ts':ts_up}
        
    elif comp=='barystatic':
        ds=ds.sel(names=['IMB_WGP'])
        ds = ds.where((ds.lat>-66) & (ds.lat<66),np.nan)

        trend = np.array(ds.SLA[0,0,:,:])
        unc = np.array(ds.SLA_UNC[0,0,:,:])
        da=da['SLF_ASL'].sel(reconstruction='ENS')
        da = da.where((da.lat>-66) & (da.lat<66),np.nan)
        ts = np.array(da.data)
    
    dic[comp] = {'trend':trend,
                'unc':unc,
                'ts':ts
                }

lat=np.array(ds.lat)
lon=np.array(ds.lon)

#%% sum of comps and residual
datasets = ['steric','barystatic','dynamic']
das_unc = []
das_trend = []
das_ts = []
for key in datasets:
    das_unc.append(dic[key]['unc'])
    das_trend.append(dic[key]['trend'])
    das_ts.append(dic[key]['ts'])
    
sum_comps_trend = sum_linear(das_trend)
sum_comps_unc = sum_linear(das_unc)
sum_comps_ts = sum_linear(das_ts)

res_trend = sum_linear([dic['alt']['trend'], sum_comps_trend], how='subtract')
res_unc = sum_square([dic['alt']['unc'], sum_comps_unc], how='subtract')
res_ts = sum_linear([dic['alt']['ts'], sum_comps_ts], how='subtract')

dic['res']={'trend':res_trend,'unc':res_unc,'ts':res_ts}
dic['sum']={'trend':sum_comps_trend,'unc':sum_comps_unc,'ts':sum_comps_ts}

#% % make list with datasets
datasets = ['alt','steric','barystatic','dynamic','sum','res']
titles = [r"$\eta_{sat(cor)}$",
          r"$\eta_{SSL}$",
          r"$\eta_{BSL}$",
          r"$\eta_{DSL}$",
          r"$\sum(\eta_{SSL}+\eta_{BSL}+\eta_{DSL})$", 
          r"$\eta_{sat} -  \eta_{\sum}$"]
# \eta_{obs} = \eta_{SSL} = \eta_{BSL} + \eta_{DSL} 
# plt.title(r"$\eta$")
das_unc = []
das_trend = []
das_ts = []
for key in datasets:
    das_unc.append(dic[key]['unc'])
    das_trend.append(dic[key]['trend'])
    das_ts.append(dic[key]['ts'])
    

#%% plot trends for each component
clim=5
plot_map_subplots( das_trend,
             plot_type = 'pcolor',
             lon=lon,lat=lat,
             cmap=cmap_trend,
             cmin=-clim,cmax=clim,
             titles=titles,
             clabel='Trend \nmm/yr',
             lon0=210, offset_y = -0.15,
             fontsize=25,
             nrow=3,ncol=2
             )

#%% lower panel:
fontsize=25
plt.figure(figsize=(15,10),dpi=300)
ax2 = plt.subplot(111)
for idata,data in enumerate(das_ts):
    data = data*landmask
    mu = np.nanmean(data,axis=(1,2))
    out = sl.get_ts_trend(tdec,mu,plot=False)
    tr = np.round(out[0],2)
    if tr==0:
        tr=0.00
    ax2.plot(tdec, mu - np.nanmean(mu[144:276]),
            label='{}: {:.2f} mm/yr'.format(titles[idata],tr),
            linewidth=3)
plt.title('Global Mean Sea Level',fontsize=fontsize)
plt.ylabel('mm',fontsize=fontsize-5)
plt.xlabel('time',fontsize=fontsize-5)

#. plt.legend(fontsize=fontsize-5)
ax2.legend(loc='lower center', bbox_to_anchor=(0.5,-0.2),
          ncol=3, 
           fancybox=True, 
           shadow=True,
           fontsize=fontsize-8)


#%% plot unc
clim=3
plot_map_subplots( das_unc,
             plot_type = 'pcolor',
             lon=lon,lat=lat,
             cmap=cmap_unc,
             cmin=0,cmax=clim,
             titles=titles,
             clabel='Uncertainty \nmm/yr',
             lon0=210, offset_y = -0.15,
             fontsize=25,
             nrow=3,ncol=2
             )
#%% plot only trend and residuals > unc
das_sig = []
for key in datasets:
    tr = np.array(dic[key]['trend'])
    unc = np.array(dic[key]['unc'])
    tr[np.abs(tr)<unc] = np.nan
    das_sig.append(tr)
clim=5
plot_map_subplots( das_sig,
             plot_type = 'pcolor',
             lon=lon,lat=lat,
             cmap=cmap_trend,
             cmin=-clim,cmax=clim,
             titles=titles,
             clabel='Signifcant Trend (larger than uncertainty) \nmm/yr',
             lon0=210, offset_y = -0.15,
             fontsize=25,
             nrow=3,ncol=2
             )

res_sg = das_sig[-1]
res_tr = das_trend[-1]
n_cells = len(res_tr[np.isfinite(res_tr)])
n_open = len(res_sg[np.isfinite(res_sg)])
n_closed = n_cells - n_open
perc_closed = n_closed / n_cells * 100
perc_open = n_open / n_cells * 100

print('Out of {}, {} have a residual lower than the uncertainty'.format(n_cells, n_closed))
print('That is, {} % of the budget is closed'.format(perc_closed))
print('and {} of the budget is unsolved ({} cells)'.format(perc_open,n_open))
#%% unc = 0.3
das_sig = []
for key in datasets:
    tr = np.array(dic[key]['trend'])
    unc = np.array(dic[key]['unc'])
    unc=0.3
    tr[np.abs(tr)<unc] = np.nan
    das_sig.append(tr)
clim=5
plot_map_subplots( das_sig,
             plot_type = 'pcolor',
             lon=lon,lat=lat,
             cmap=cmap_trend,
             cmin=-clim,cmax=clim,
             titles=titles,
             clabel='Signifcant Trend (larger than uncertainty) \nmm/yr',
             lon0=210, offset_y = -0.15,
             fontsize=25,
             nrow=3,ncol=2
             )
res_sg = das_sig[-1]
res_tr = das_trend[-1]
n_cells = len(res_tr[np.isfinite(res_tr)])
n_open = len(res_sg[np.isfinite(res_sg)])
n_closed = n_cells - n_open
perc_closed = n_closed / n_cells * 100
perc_open = n_open / n_cells * 100

print('Out of {}, {} have a residual lower than the uncertainty'.format(n_cells, n_closed))
print('That is, {} % of the budget is closed'.format(perc_closed))
print('and {} of the budget is unsolved ({} cells)'.format(perc_open,n_open))
#%%
# from matplotlib import colors
# #%%
# c0=0
# cmax_unc=3
# plt.figure()
# ax2=plt.subplot(111)
# d  = das_unc[-1]
# # plt.hist(das_unc[-1])
# d=np.hstack(d)
# # equal_area = [np.repeat(d[i],np.round(area[i]/1)) for i in range(1,len(d)) ]
# # equal_area = np.concatenate( equal_area, axis=0 )
# # d=equal_area
# bins= np.arange(np.nanmin(d)-1,np.nanmax(d)+1,0.1)

# n, bins, patches = plt.hist(d, bins=bins, # bins=int(rg/bw), 
#                             density=True, 
#                             facecolor='#2ab0ff', edgecolor='#e0e0e0', 
#                             linewidth=0.5, alpha=0.5)

# #norm = colors.Normalize(0,1)
# norm = colors.Normalize(c0,cmax_unc)

# for bin,patch in zip(bins,patches):
#     color=cmap_unc(norm(bin))
#     patch.set_facecolor(color)
# ax2.set_xlim(c0,cmax_unc)
# #%%

# df= pd.DataFrame( {'unc':np.hstack(das_unc[-1]),
#                    'trend':np.hstack(das_trend[-1])})
# N_cells = len(df)
# df.dropna(inplace=True)
# N_cells_oceans = len(df)
# df.reset_index(inplace=True)


# #% %

# data_x = df.index
# data_y = df.unc
# data_y_sc = np.abs(df.trend)
# N_cells_closed = len(data_y_sc[data_y_sc>data_y])
# print('Out of {}, {} have a residual lower than the uncertainty'.format(N_cells_oceans, N_cells_closed))
# print('That is, {} % '.format(N_cells_closed * 100 /N_cells_oceans ))
# #%%
# data_x = df.index[0:6500]
# data_y = df.unc[0:6500]
# data_y_sc = df.trend[0:6500]

# res_closed = len(data_y_sc[data_y_sc>data_y])
# res_max = len(data_y_sc)
# print()


# fig, ax = plt.subplots(figsize=(25, 4))
# # plt.xticks(data_x)    
# plt.ylabel("mm/yr")
# plt.title('Budget Residuals: {} out of {} closed'.format(res_closed,res_max))

# ax.set_xlim(np.nanmin(data_x),np.nanmax(data_x))
# # bars
# my_cmap = cmap_unc
# data_y_normalized = [x / max(data_y) for x in data_y]
# colors = my_cmap(data_y_normalized)
# sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(0,3))
# sm.set_array([])
# cbar_ax = fig.add_axes([0.96, 0.15, 0.01, 0.7])
# cbar = plt.colorbar(sm,cax=cbar_ax)
# cbar.set_label('Uncertainty', rotation=270,labelpad=15)
# rects = ax.bar(data_x, data_y, 
#                 # color=colors, 
#                 alpha=0.5,
#                edgecolor=None)
# #% %
# # scatter
# my_cmap = cmap_trend
# data_y_normalized = [x / max(data_y_sc) for x in data_y]
# colors = my_cmap(data_y_sc)
# sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(-5,5)
#                     )
# sm.set_array([])
# cbar_ax2 = fig.add_axes([0.91, 0.15, 0.01, 0.7])
# cbar2 = plt.colorbar(sm,cax=cbar_ax2)
# cbar2.set_label('Trend', rotation=270,labelpad=15)
# #% %
# sc = ax.scatter(data_x,np.abs(data_y_sc),
#                 s=5,
#                 c=cmap_trend(data_y_sc))


# # plt.savefig("bar_chart_with_colorbar_03.png", bbox_inches='tight')

# plt.show()

#%% Clusters

#%% SOM 19 regions
path = '//Volumes/LaCie_NIOZ/budget/regions/som/'
file= 'som_3x3_alt_1993_2019_n10_sig2_ep_atlantic_indopacific'
ds=xr.open_dataset(path+file+'.nc')

mask_clusters = np.array(ds.mask)
som = {'mask':mask_clusters}
n_clusters = len(np.unique(mask_clusters[np.isfinite(mask_clusters)]))
plot_map(mask_clusters,cmax=n_clusters,lon0=210,title='SOM Clusters',clabel='cluster number')

# mask_comb = np.array(mask_clusters)
x=4
y=5
plt.figure(figsize=(20,10))
for i in range(0,int(np.nanmax(mask_clusters))):
    icluster = i+1
    ax = plt.subplot(y,x,icluster, projection = ccrs.Robinson(central_longitude=210))
    ax.set_global()
    mask=np.array(mask_clusters)
    mask[np.where(mask!=icluster)]=np.nan
    mm = ax.pcolormesh(ds.lon,\
                       ds.lat,\
                       mask,
                       vmin=0, vmax=x*y, 
                       transform=ccrs.PlateCarree(),
                       #cmap='Spectral_r'
                       cmap='jet'
                      )
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', 
                                                edgecolor='gray', facecolor='papayawhip'))
    plt.title('Cluster {}'.format(icluster))
plt.show()

som['n']=n_clusters

#% % SOM 19 residuals
mat = np.zeros((n_clusters,len(datasets)))
mat2 = np.zeros((n_clusters,len(datasets)))
mat3 = np.zeros((n_clusters))

n = []

for j,label in enumerate(datasets):
    tmp = np.zeros((n_clusters,180,360))
    tmp2 = np.zeros((n_clusters,180,360))
    tmp3 = np.zeros((n_clusters,180,360))
    
    test = np.zeros((n_clusters,180,360))
    for i in range(n_clusters):
        icluster = i+1
        mask=np.array(mask_clusters)
        mask[np.where(mask!=icluster)]=np.nan
        mask[np.isfinite(mask)]=1
        
        # tmp[i,mask==1] = cluster_mean(np.array(dic[label]['trend']),mask, lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
        # tmp2[i,mask==1] = cluster_mean(np.array(dic[label]['unc']),mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
        test[i,mask==1] = icluster
        mat[i,j] = cluster_mean(np.array(dic[label]['trend']),mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
        mat2[i,j] = cluster_mean(np.array(dic[label]['unc']),mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
        tmp[i,mask==1] = mat[i,j]
        tmp2[i,mask==1] = mat2[i,j]
        if label =='res':
               mat3[i] = cluster_mean(np.array( np.abs(dic['alt']['unc'] -dic['sum']['unc']) ),mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
               tmp3[i,mask==1]=mat3[i]
               
    tr = np.sum(tmp,axis=0)
    tr[tr==0] = np.nan
    unc = np.sum(tmp2,axis=0)
    unc[unc==0] = np.nan
    
    som[label] = {'trend':tr,
                 'unc':unc,
                 }
    if label=='res':
        unc2 = np.sum(tmp3,axis=0)
        unc2[unc2==0] = np.nan
        som[label]['unc2']=unc2
        
som['res']['unc2_cl']=mat3


df_som = pd.DataFrame ({'cluster_n': np.unique(mask_clusters[np.isfinite(mask_clusters)]) })
for j,label in enumerate(datasets):
    df_som['{}_tr'.format(label)] = mat[:,j]
    df_som['{}_unc'.format(label)] = mat2[:,j]
df_som   

som['df']=df_som
dic['som']=som

#% %
#%% dmaps regions

ds=xr.open_dataset('/Volumes/LaCie_NIOZ/budget/regions/dmaps/dmaps_k5_k23.nc')

mask_clusters = np.array(ds.mask[0,:,:])
dmap = {'mask':mask_clusters}
n_clusters = len(np.unique(mask_clusters[np.isfinite(mask_clusters)]))
plot_map(mask_clusters,cmax=n_clusters,lon0=210,cmap='prism',title='dMAPS k5 Clusters',clabel='cluster number')

# mask_comb = np.array(mask_clusters)
# x=5
# y=20
# plt.figure(figsize=(15,35))
# for i in range(0,int(np.nanmax(mask_clusters))):
#     icluster = i+1
#     ax = plt.subplot(y,x,icluster, projection = ccrs.Robinson(central_longitude=210))
#     ax.set_global()
#     mask=np.array(mask_clusters)
#     mask[np.where(mask!=icluster)]=np.nan
#     mm = ax.pcolormesh(ds.lon,\
#                        ds.lat,\
#                        mask,
#                        vmin=0, vmax=x*y, 
#                        transform=ccrs.PlateCarree(),
#                        #cmap='Spectral_r'
#                        cmap='jet'
#                       )
#     ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', 
#                                                 edgecolor='gray', facecolor='papayawhip'))
#     plt.title('Cluster {}'.format(icluster))
# plt.show()


#% % dmaps residuals
mat = np.zeros((n_clusters,len(datasets)))
mat2 = np.zeros((n_clusters,len(datasets)))
mat3 = np.zeros((n_clusters))
n = []
dmap['n']=n_clusters

for j,label in enumerate(datasets):
    tmp = np.zeros((n_clusters,180,360))
    tmp2 = np.zeros((n_clusters,180,360))
    tmp3 = np.zeros((n_clusters,180,360))
    
    test = np.zeros((n_clusters,180,360))
    for i in range(n_clusters):
        icluster = i+1
        mask=np.array(mask_clusters)
        mask[np.where(mask!=icluster)]=np.nan
        mask[np.isfinite(mask)]=1
        
        # tmp[i,mask==1] = cluster_mean(np.array(dic[label]['trend']),mask, lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
        # tmp2[i,mask==1] = cluster_mean(np.array(dic[label]['unc']),mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
        test[i,mask==1] = icluster
        mat[i,j] = cluster_mean(np.array(dic[label]['trend']),mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
        mat2[i,j] = cluster_mean(np.array(dic[label]['unc']),mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
        tmp[i,mask==1] = mat[i,j]
        tmp2[i,mask==1] = mat2[i,j]
        if label =='res':
                mat3[i] = cluster_mean(np.array( np.abs(dic['alt']['unc'] -dic['sum']['unc']) ),mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
                tmp3[i,mask==1] = mat3[i]
                
    tr = np.sum(tmp,axis=0)
    tr[tr==0] = np.nan
    unc = np.sum(tmp2,axis=0)
    unc[unc==0] = np.nan
    
    dmap[label] = {'trend':tr,
                 'unc':unc,
                 }
    if label=='res':
        unc2 = np.sum(tmp3,axis=0)
        unc2[unc2==0] = np.nan
        dmap[label]['unc2']=unc2
        
dmap['res']['unc2_cl']=mat3

df_dmap = pd.DataFrame ({'cluster_n': np.unique(mask_clusters[np.isfinite(mask_clusters)]) })
for j,label in enumerate(datasets):
    df_dmap['{}_tr'.format(label)] = mat[:,j]
    df_dmap['{}_unc'.format(label)] = mat2[:,j]
df_dmap   

dmap['df']=df_dmap
dic['dmap']=dmap

#%% Plot residual clusters
#% % make list with datasets
das_clusters = [dic['som']['res']['trend'],dic['dmap']['res']['trend'],
                # dic['som']['res']['unc'],dic['dmap']['res']['unc'],
                ]
tr = np.array(dic['som']['res']['trend'])
unc = np.array(dic['som']['res']['unc2'])
tr[np.abs(tr)<unc] = np.nan
das_clusters.append(tr)
tr = np.array(dic['dmap']['res']['trend'])
unc = np.array(dic['dmap']['res']['unc2'])
tr[np.abs(tr)<unc] = np.nan
das_clusters.append(tr)

titles = [r"SOM residuals ",r"dMAPS residuals ", 
          # r"unc", r"unc", 
          'Significant residual', "Significant residual"
                    
          ]

#% % plot trends for each component
clim=5
plot_map_subplots( das_clusters,
             plot_type = 'pcolor',
             lon=lon,lat=lat,
             cmap=cmap_trend,
             cmin=-clim,cmax=clim,
             titles=titles,
             clabel='nmm/yr',
             lon0=210, offset_y = +0.2,
             fontsize=25,
             fsize=(10,10),
             nrow=3,ncol=2
             )

#%% plot residual with relation to 0.3
das_clusters = [dic['som']['res']['trend'],dic['dmap']['res']['trend'],
                # dic['som']['res']['unc'],dic['dmap']['res']['unc'],
                ]
tr = np.array(dic['som']['res']['trend'])
unc=0.3
# unc = np.array(dic['som']['res']['unc'])

tr[np.abs(tr)<unc] = np.nan
das_clusters.append(tr)
tr = np.array(dic['dmap']['res']['trend'])
# unc = np.array(dic['dmap']['res']['unc'])
tr[np.abs(tr)<unc] = np.nan
das_clusters.append(tr)

titles = [r"SOM residuals ",r"dMAPS residuals ", 
          # r"unc", r"unc", 
          'Significant residual (0.3)', "Significant residual (0.3)"
                    
          ]

#% % plot trends for each component
clim=5
plot_map_subplots( das_clusters,
             plot_type = 'pcolor',
             lon=lon,lat=lat,
             cmap=cmap_trend,
             cmin=-clim,cmax=clim,
             titles=titles,
             clabel='nmm/yr',
             lon0=210, offset_y = +0.2,
             fontsize=25,
             fsize=(10,10),
             nrow=3,ncol=2
             )
#%%

# print()

fig = plt.figure(figsize=(15,10))
ax = plt.subplot(211)
# plt.xticks(data_x)    
key = 'som'
df= pd.DataFrame( {
                    # 'unc':np.hstack(dic['som']['df']['res_unc']),
                   'unc':dic[key]['res']['unc2_cl'],
                    'trend':np.hstack(dic[key]['df']['res_tr']) })
# N_cells = len(df)
df.dropna(inplace=True)
# N_cells_oceans = len(df)
df.reset_index(inplace=True)

data_x = df.index
data_y = df.unc
data_y_sc = np.abs(df.trend)
# N_cells_closed = len(data_y_sc[data_y_sc>data_y])
# print('Out of {}, {} have a residual lower than the uncertainty'.format(N_cells_oceans, N_cells_closed))
# print('That is, {} % '.format(N_cells_closed * 100 /N_cells_oceans ))
# #%%
# data_x = df.index[0:6500]
# data_y = df.unc[0:6500]
# data_y_sc = df.trend[0:6500]

res_closed = len(data_y_sc[np.abs(data_y_sc)<data_y])
res_max = len(data_y_sc)
plt.ylabel("mm/yr")
plt.title('SOM Residuals: {} out of {} closed'.format(res_closed,res_max))

# ax.set_xlim(np.nanmin(data_x),np.nanmax(data_x))
# bars
my_cmap = cmap_unc
data_y_normalized = [x / max(data_y) for x in data_y]
colors = my_cmap(data_y_normalized)
rects = ax.bar(data_x, data_y, 
                color=colors, 
                alpha=0.5,
                edgecolor=None)
#% %
# scatter
#% %
sc = ax.scatter(data_x,np.abs(data_y_sc),
                s=50,
                c=cmap_trend(data_y_sc))


ax = plt.subplot(212)
# plt.xticks(data_x)    
key = 'dmap'
df= pd.DataFrame( {
                    # 'unc':np.hstack(dic['som']['df']['res_unc']),
                   'unc':dic[key]['res']['unc2_cl'],
                    'trend':np.hstack(dic[key]['df']['res_tr']) })
# N_cells = len(df)
df.dropna(inplace=True)
# N_cells_oceans = len(df)
df.reset_index(inplace=True)

data_x = df.index
data_y = df.unc
data_y_sc = np.abs(df.trend)
# N_cells_closed = len(data_y_sc[data_y_sc>data_y])
# print('Out of {}, {} have a residual lower than the uncertainty'.format(N_cells_oceans, N_cells_closed))
# print('That is, {} % '.format(N_cells_closed * 100 /N_cells_oceans ))
# #%%
# data_x = df.index[0:6500]
# data_y = df.unc[0:6500]
# data_y_sc = df.trend[0:6500]

res_closed = len(data_y_sc[np.abs(data_y_sc)<data_y])
res_max = len(data_y_sc)
plt.ylabel("mm/yr")
plt.title('dMAPS Residuals: {} out of {} closed'.format(res_closed,res_max))

# ax.set_xlim(np.nanmin(data_x),np.nanmax(data_x))
# bars
my_cmap = cmap_unc
data_y_normalized = [x / max(data_y) for x in data_y]
colors = my_cmap(data_y_normalized)
rects = ax.bar(data_x, data_y, 
                color=colors, 
                alpha=0.5,
                edgecolor=None)
#% %
# scatter
#% %
sc = ax.scatter(data_x,np.abs(data_y_sc),
                s=50,
                c=cmap_trend(data_y_sc))


plt.xlabel('cluster number')

#% color bars
sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(0,3))
sm.set_array([])
cbar_ax = fig.add_axes([0.98, 0.15, 0.02, 0.7])
cbar = plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Uncertainty', rotation=270,labelpad=15)


my_cmap = cmap_trend
colors = my_cmap(data_y_sc)
sc = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(-5,5)
                    )
sc.set_array([])
cbar_ax2 = fig.add_axes([0.91, 0.15, 0.02, 0.7])
cbar2 = plt.colorbar(sc,cax=cbar_ax2)
cbar2.set_label('Trend', rotation=270,labelpad=15)



# plt.savefig("bar_chart_with_colorbar_03.png", bbox_inches='tight')

plt.show()

#%% SOM 4x4 regions

path = '//Volumes/LaCie_NIOZ/budget/regions/som/'
file= 'som_4x4_alt_1993_2019_n10_sig2_ep_world'
ds=xr.open_dataset(path+file+'.nc')

mask_clusters = np.array(ds.regions)
som = {'mask':mask_clusters}
n_clusters = len(np.unique(mask_clusters[np.isfinite(mask_clusters)]))

plot_map(mask_clusters,cmax=n_clusters,lon0=210,title='SOM Clusters',clabel='cluster number')

# mask_comb = np.array(mask_clusters)
x=4
y=5
plt.figure(figsize=(20,10))
for i in range(0,int(np.nanmax(mask_clusters))):
    icluster = i+1
    ax = plt.subplot(y,x,icluster, projection = ccrs.Robinson(central_longitude=210))
    ax.set_global()
    mask=np.array(mask_clusters)
    mask[np.where(mask!=icluster)]=np.nan
    mm = ax.pcolormesh(ds.lon,\
                       ds.lat,\
                       mask,
                       vmin=0, vmax=x*y, 
                       transform=ccrs.PlateCarree(),
                       #cmap='Spectral_r'
                       cmap='jet'
                      )
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', 
                                                edgecolor='gray', facecolor='papayawhip'))
    plt.title('Cluster {}'.format(icluster))
plt.show()
# cluster 10 is empty:
n_clusters = len(np.arange(0,int(np.nanmax(mask_clusters))))
som['n']=n_clusters
#%% SOM 19 residuals
mat = np.zeros((n_clusters,len(datasets)))
mat2 = np.zeros((n_clusters,len(datasets)))
mat = np.zeros((n_clusters,len(datasets)))
mat2 = np.zeros((n_clusters,len(datasets)))
mat3 = np.zeros((n_clusters))

n = []

for j,label in enumerate(datasets):
    tmp = np.zeros((n_clusters,180,360))
    tmp2 = np.zeros((n_clusters,180,360))
    tmp3 = np.zeros((n_clusters,180,360))
    
    test = np.zeros((n_clusters,180,360))
    for i in range(n_clusters):
        icluster = i+1
        mask=np.array(mask_clusters)
        mask[np.where(mask!=icluster)]=np.nan
        mask[np.isfinite(mask)]=1
        
        # tmp[i,mask==1] = cluster_mean(np.array(dic[label]['trend']),mask, lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
        # tmp2[i,mask==1] = cluster_mean(np.array(dic[label]['unc']),mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
        test[i,mask==1] = icluster
        mat[i,j] = cluster_mean(np.array(dic[label]['trend']),mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
        mat2[i,j] = cluster_mean(np.array(dic[label]['unc']),mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
        tmp[i,mask==1] = mat[i,j]
        tmp2[i,mask==1] = mat2[i,j]
        if label =='res':
               mat3[i] = cluster_mean(np.array( np.abs(dic['alt']['unc'] -dic['sum']['unc']) ),mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
               tmp3[i,mask==1]=mat3[i]
               
    tr = np.sum(tmp,axis=0)
    tr[tr==0] = np.nan
    unc = np.sum(tmp2,axis=0)
    unc[unc==0] = np.nan
    
    som[label] = {'trend':tr,
                 'unc':unc,
                 }
    if label=='res':
        unc2 = np.sum(tmp3,axis=0)
        unc2[unc2==0] = np.nan
        som[label]['unc2']=unc2
        
som['res']['unc2_cl']=mat3
    

df_som = pd.DataFrame ({'cluster_n': np.arange(0,int(np.nanmax(mask_clusters))) })
for j,label in enumerate(datasets):
    df_som['{}_tr'.format(label)] = mat[:,j]
    df_som['{}_unc'.format(label)] = mat2[:,j]
df_som   

som['df']=df_som
dic['som_4x4']=som

#% % make list with datasets
das_clusters = [dic['som_4x4']['res']['trend'],# dic['dmap']['res']['trend'],
                # dic['som']['res']['unc'],dic['dmap']['res']['unc'],
                ]
tr = np.array(dic['som_4x4']['res']['trend'])
unc = np.array(dic['som_4x4']['res']['unc'])
tr[np.abs(tr)<unc] = np.nan
das_clusters.append(tr)
# tr = np.array(dic['dmap']['res']['trend'])
# unc = np.array(dic['dmap']['res']['unc'])
# tr[np.abs(tr)<unc] = np.nan
# das_clusters.append(tr)

titles = [r"SOM residuals ",# r"dMAPS residuals ", 
          # r"unc", r"unc", 
          'Significant residual', # "Significant residual"
                    
          ]

#% % plot trends for each component
clim=5
plot_map_subplots( das_clusters,
             plot_type = 'pcolor',
             lon=lon,lat=lat,
             cmap=cmap_trend,
             cmin=-clim,cmax=clim,
             titles=titles,
             clabel='nmm/yr',
             lon0=210, offset_y = -0.2,
             fontsize=25,
             fsize=(15,10),
             nrow=2,ncol=1
             )

#%% Budget sensitivity
#%% by budget component
combs = [['steric','barystatic','dynamic'],
            ['steric','barystatic'],
            ['steric_up','barystatic']
                ]
budgets=['1_deg','som','dmap']
combos=[]
matrix = np.zeros((len(combs),len(budgets)))
for ic,dataset in enumerate(combs):
    if len (dataset)==2:
        comb = '{}{}'.format(dataset[0],dataset[1])
    elif len(dataset)==3:
        comb = '{}{}{}'.format(dataset[0],dataset[1],dataset[2])
    combos.append(comb)
    das_unc = []
    das_trend = []
    das_ts = []
    for key in dataset:
        das_unc.append(dic[key]['unc'])
        das_trend.append(dic[key]['trend'])
        das_ts.append(dic[key]['ts'])
        
    sum_comps_trend = sum_linear(das_trend)
    sum_comps_unc = sum_linear(das_unc)
    sum_comps_ts = sum_linear(das_ts)
    
    res_trend = sum_linear([dic['alt']['trend'], sum_comps_trend], how='subtract')
    res_unc = sum_square([dic['alt']['unc'], sum_comps_unc], how='subtract')
    res_ts = sum_linear([dic['alt']['ts'], sum_comps_ts], how='subtract')

    dic[comb] = {'res':{'trend':res_trend,'unc':res_unc,'ts':res_ts},
                  'sum':{'trend':sum_comps_trend,'unc':sum_comps_unc,'ts':sum_comps_ts}
                  }
    for ib, budget in enumerate(budgets):
        if budget=='1_deg':
            n_cells = len(res_trend[np.isfinite(res_tr)])
            res_sig = np.array(res_trend)
            res_sig[res_sig<res_unc] = np.nan
            n_open = len(res_sig[np.isfinite(res_sig)])
            n_close = n_cells-n_open
        else:
            # here
            mask_clusters = dic[budget]['mask']
            n_clusters = dic[budget]['n']
            
            #% %cluster residuals
            mat = np.zeros((n_clusters))
            mat2 = np.zeros((n_clusters,))
            mat3 = np.zeros((n_clusters))
            
            n = []
            
            tmp = np.zeros((n_clusters,180,360))
            tmp2 = np.zeros((n_clusters,180,360))
            tmp3 = np.zeros((n_clusters,180,360))
                
            for i in range(n_clusters):
                icluster = i+1
                mask=np.array(mask_clusters)
                mask[np.where(mask!=icluster)]=np.nan
                mask[np.isfinite(mask)]=1
                
                # tmp[i,mask==1] = cluster_mean(np.array(dic[label]['trend']),mask, lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
                # tmp2[i,mask==1] = cluster_mean(np.array(dic[label]['unc']),mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
                mat[i] = cluster_mean(res_trend,mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
                mat2[i] = cluster_mean(res_unc,mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
                tmp[i,mask==1] = mat[i]
                tmp2[i,mask==1] = mat2[i]

                mat3[i] = cluster_mean(np.array( np.abs(dic['alt']['unc'] - sum_comps_unc) ),mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
                tmp3[i,mask==1]=mat3[i]
                           
            tr = np.sum(tmp,axis=0)
            tr[tr==0] = np.nan
            unc = np.sum(tmp2,axis=0)
            unc[unc==0] = np.nan
                
            dic[comb][budget] = {'trend':tr,
                          'unc':unc,
                          }
            unc2 = np.sum(tmp3,axis=0)
            unc2[unc2==0] = np.nan
            dic[comb][budget]['unc2']=unc2
                    
            dic[comb][budget]['unc2_cl']=mat3
                
            
            n_cells = n_clusters
            res_sig = np.array(mat)
            res_sig[res_sig<mat2] = np.nan
            n_close = len(res_sig[np.isnan(res_sig)])
        matrix[ic,ib] = (n_close/n_cells) * 100
        # print(comb)
        # print(budgets[ib])
        # print(matrix[ic,ib])
        # print('')
#%%
df= pd.DataFrame(matrix,columns=['1deg','som','dmap'], 
                 index=['steric, barystatic, dynamic', 'steric, barystatic', 'steric up, barystatic'])
# df['combination']=combos
print(df)
#%%
fig = plt.figure(dpi=300)
ax = plt.subplot(111)
df.plot.bar(rot=15,ax=ax).legend(loc='lower right')
plt.ylabel('% closed budget')
plt.show()
#%%
fig = plt.figure()
ax = plt.subplot(111)
df.plot.barh(rot=45,ax=ax).legend(loc='upper right')
ax.set_xlabel('% closed budget')
plt.show()

#%% sensitivity to budget combination 
da = xr.open_dataset('/Volumes/LaCie_NIOZ/data/budget/combinations.nc')
budget = 'som'
mask_clusters = dic[budget]['mask']
n_clusters = dic[budget]['n']

n_pos = len(da.comb)
cluster_combos_res = np.zeros((n_clusters,n_pos))
cluster_combos_unc = np.zeros((n_clusters,n_pos))
cluster_combos_sig = np.zeros((n_clusters,n_pos))

for ipos in range(n_pos):
    res=np.array(da.res[ipos])
    unc = np.array(da.unc[ipos])
    mat = np.zeros((n_clusters))
    mat2 = np.zeros((n_clusters,))
    # mat3 = np.zeros((n_clusters))
    
    for i in range(n_clusters):
        icluster = i+1
        mask=np.array(mask_clusters)
        mask[np.where(mask!=icluster)]=np.nan
        mask[np.isfinite(mask)]=1
        
        # tmp[i,mask==1] = cluster_mean(np.array(dic[label]['trend']),mask, lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
        # tmp2[i,mask==1] = cluster_mean(np.array(dic[label]['unc']),mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
        mat[i] = cluster_mean(res,mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
        mat2[i] = cluster_mean(unc,mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
        # mat3[i] = cluster_mean(np.array( np.abs(dic['alt']['unc'] - sum_comps_unc) ),mask,lat=np.array(ds.lat),lon=np.array(ds.lon),norm=False )
        cluster_combos_res[i,ipos] = mat[i]
        cluster_combos_unc[i,ipos] = mat2[i]
        if mat[i]<mat2[i]:
            cluster_combos_sig[i,ipos] = 0 # zero if closed
        else:
            cluster_combos_sig[i,ipos] = 1 # 1 if open
#%% plot
plt.figure(dpi=300)
ax=plt.subplot(111)
alpha=0.1
x=range(n_clusters)
for ipos in range(n_pos):
    if np.any(np.isfinite(cluster_combos_res[:,ipos])):
        col = ['salmon' if sig==1 else 'mediumslateblue' for sig in cluster_combos_sig[:,ipos]]
        sc=plt.scatter(x,cluster_combos_res[:,ipos],
                        c=col,
                        # c=cluster_combos_sig[:,ipos],
                       alpha=alpha/2)

# plot just one to get lengend
ind=np.where(cluster_combos_sig[:,ipos]==1)[0][0]
plt.scatter(x[ind],cluster_combos_res[ind,ipos],
            c=col[ind],
            alpha=alpha,label='open')
ind=np.where(cluster_combos_sig[:,ipos]==0)[0][0]
plt.scatter(x[ind],cluster_combos_res[ind,ipos],
            c=col[ind],
            alpha=alpha,label='closed')

## ENS
names = np.array(da.names)
ipos =[]
ipos = [i for i in range(len(da.comb)) if np.all(names[i,:]==['ENS'])][0]
plt.scatter(x,cluster_combos_res[:,ipos],
            # c='black',
            marker='s',
            c = ['salmon' if sig==1 else 'mediumslateblue' for sig in cluster_combos_sig[:,ipos]],
            # label='ENS'
            )
ind=np.where(cluster_combos_sig[:,ipos]==1)[0][0]
plt.scatter(x[ind],cluster_combos_res[ind,ipos],
            marker='s',
            c=col[ind],
            # alpha=alpha,
            label='ENS - open')
ind=np.where(cluster_combos_sig[:,ipos]==0)[0][0]
plt.scatter(x[ind],cluster_combos_res[ind,ipos],
            c=col[ind],
            marker='s',
            # alpha=alpha,
            label='ENS - closed')
plt.legend(ncol=2)
plt.ylabel('Residual trend mm/yr')
plt.xlabel('SOM cluster number')

#%%


