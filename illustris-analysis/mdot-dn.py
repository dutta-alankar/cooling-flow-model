# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 12:04:17 2021

@author: alankar
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42
#matplotlib.use('pdf')
from decimal import Decimal

def fexp(number):
    (sign, digits, exponent) = Decimal(number).as_tuple()
    return len(digits) + exponent - 1

def fman(number):
    return Decimal(number).scaleb(-fexp(number)).normalize()

Msun = 2e33
yr = 365*24*60**2
mp = 1.6726219e-24
kB = 1.380649e-16
pc = 3.086e18
kpc = 1e3*pc

X, Y, Z = 0.7154, 0.2703, 0.0143
Z *= 10**-0.2
X = 1.-(Y+Z)
mu = 1/(2*X+0.75*Y+0.5625*Z)
mue = 2./(1+X)
mui = 1/(1/mu-1/mue)
muH = 1/X

gamma = 5./3
gamma_m = 1.03
sizes = [0.5, 1.0, 1.5, 2.0] #kpc

K = 4*np.pi
q = 2

mdot_50, mdot_16, mdot_84 = [], [], []
Be, dn_temperature = [], []
dn_radius = None
for select in range(len(sizes)):
    r_crit = sizes[select]
    start = select*3+1
    print('Choosing cloud of size %.1f kpc'%sizes[select])
    
    #----------------------------------
    dn_data = np.loadtxt('../dylan-nelson/fig10_TNG50-1_z0.5_h8_nH.txt', skiprows=5)
    dn_radius = dn_data[:,0]
    #print(dn_radius)
    dn_nH_50 = 10**dn_data[:,start]
    dn_nH_16, dn_nH_84 = dn_data[:,start+1],dn_data[:,start+2]
    
    #----------------------------------
    dn_data = np.loadtxt('../dylan-nelson/fig10_TNG50-1_z0.5_h8_vel_rel.txt', skiprows=5)
    dn_radius = dn_data[:,0]
    dn_vel_50 = dn_data[:,start]
    dn_vel_16, dn_vel_84 = dn_data[:,start+1],dn_data[:,start+2]
    
    mdot_16.append( -K*(dn_radius*kpc)**q*(dn_nH_16*mu*mp)*(dn_vel_16*1e5)/(Msun/yr) )
    mdot_50.append( -K*(dn_radius*kpc)**q*(dn_nH_50*mu*mp)*(dn_vel_50*1e5)/(Msun/yr) )
    mdot_84.append( -K*(dn_radius*kpc)**q*(dn_nH_84*mu*mp)*(dn_vel_84*1e5)/(Msun/yr) )


sizes = np.hstack(([0],sizes))
plt.figure(figsize=(13,10))
for cloud in range(1,len(sizes)): 
    plt.plot(dn_radius, mdot_50[cloud-1], linewidth=5, 
             label=r'$\rm R_{cl}=%.1f \ kpc\ -\ %.1f \ kpc$'%(sizes[cloud-1], 
             sizes[cloud]) )
    plt.fill_between(dn_radius, mdot_16[cloud-1], mdot_84[cloud-1],
                    alpha=0.15 )
    
plt.yscale('log')  
plt.ylim(ymin=1e-5)  
plt.xlim(xmin=0, xmax=5.)
#plt.title('Mass flow around IllustrisTNG50-1  simulated clouds', size=28)
plt.ylabel(r'$\dot{\rm M}$ [$\rm M_\odot yr^{-1}$]',size=28)
plt.xlabel(r'distance [$\rm kpc$]',size=28)
plt.tick_params(axis='both', which='major', labelsize=24, direction="out", pad=5)
plt.tick_params(axis='both', which='minor', labelsize=24, direction="out", pad=5)
plt.grid()
#plt.ylim(ymin=-1.5e2,ymax=1.5e2)
plt.legend(loc='upper left', prop={'size': 24},framealpha=0.3, shadow=False, fancybox=True)
plt.savefig('mdot-tng.png', transparent=True, bbox_inches='tight')
plt.show()
