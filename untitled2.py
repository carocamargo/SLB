#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 11:22:00 2022

@author: ccamargo
"""

import numpy as np
import matplotlib.pyplot as plt
#%%
import datetime as dt
import xarray as xr
path = '/Volumes/LaCie_NIOZ/data/MSL/'
file = 'MSL_Serie_MERGED_Global_AVISO_GIA_Adjust_Filter2m.nc'
ds=xr.open_dataset(path+file)
rate = np.array(ds.msl*1000) #mm
time = np.array(ds.time)
def get_dectime(time):
    from datetime import datetime
    import numpy as np
    t = [datetime.utcfromtimestamp((t- np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')) 
         for t in time]
    t = [t.timetuple().tm_year + (t.timetuple().tm_yday/365) for t in t]
    return np.array(t)
t = get_dectime(time)
rate_unc = 0.3
#%%
sr = 1000
n = 300
period = n/sr
time = np.linspace(0, period, n)
signal_pure = 100*np.sin(2*np.pi*13*time)

sigma = np.sqrt(sr/2)
noise = np.random.normal(0, sigma, n)
plt.plot(signal_pure,label='pure')
plt.plot(signal_pure+noise,label='pure+noise',alpha=0.5)
plt.legend()

#%%
def create_noise(mu,sigma,t):
    noise = np.random.normal(mu,sigma,t)
    return noise
def MC_perturb(rate,N):
    t = len(rate)
    samples = np.zeros((N,t))
    for i in range(N):
        samples[i] = create_noise(0,1,t) * rate
    perturbed = np.cumsum(samples,axis=1)
    # perturbed = np.mean(samples,axis=0)
    
    return perturbed
noise_ens = MC_perturb(noise, 1000)
#%%
def MC_estimation_slope(M):
    MC_betas = []
    MC_samples = {}

    for i in range(M):
        # randomly sampling from normal distribution as error terms
        e = np.random.normal(mu, sigma, n)
        # generating independent variable by making sure the variance in X is larger than the variance in error terms
        X = 9 * np.random.normal(mu, sigma, n)
        # population distribution using the assumd parameter values alpha/beta
        Y = (alpha + beta * X + e)
        
        # running OLS regression for getting slope parameters
        model = sm.OLS(Y.reshape((-1, 1)), X.reshape((-1, 1)))
        ols_result = model.fit()
        coeff = ols_result.params
        
        MC_samples[i] = Y
        MC_betas.append(coeff)
    MC_beta_hats = np.array(MC_betas).flatten()
    return(MC_samples, MC_beta_hats)
    
MC_samples, MC_beta_hats = MC_estimation_slope(M = 10000)
beta_hat_MC = np.mean(MC_beta_hats)
#%%
def perturb_estimate(rate,rate_unc,n=1000):
    '''
    Generate an ensemble of time series by perturbing the rate with random normal noise,
    and then integrating teh rates to obtain the time series

    Parameters
    ----------
    rate : TYPE
        rate of the process
    rate_unc : TYPE
        uncertainty if the rate
    n : TYPE, optional
        Number of ensemble members. The default is 1000.

    Returns
    -------
    rate_ens:
        ensemble time series
    tseries_ens:
        integrated time series
    '''
    # Correct GRACE estimate with GIA ensemble member and perturb with measurement noise
    dt_now = dt.datetime.now().replace(microsecond=0)
    
    np.random.seed() # Reset random number generator 
    rand_normal = np.random.normal(0,1,rate.shape) # random noise
    rate_ens = np.zeros((n,np.array(rate.shape)[0]))
    for i in range(n):
        rate_ens[i] = np.array(rate + (rand_normal*rate_unc))    # Random pertubations from GRACE uncertainty
    ewh_ptb = grace['ewh'] + mscn.mascon2grid_3d(rnd_mscn, settings) - (gia['ewh']*(settings['time']-settings['time'].mean())[:,np.newaxis,np.newaxis])
    print('    ens '+str(ens+1)+' Perturbing GRACE mass:', (dt.datetime.now().replace(microsecond=0) - dt_now))
    return(ewh_ptb)