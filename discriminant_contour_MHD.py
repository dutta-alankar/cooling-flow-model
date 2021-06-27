#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 23:49:10 2020

@author: prateek
contourplot of discriminant in beta_0-Lambda_T plane
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rc('lines', linewidth=3, color='r')
plt.rcParams.update({'font.size': 16})

q=2; g=5/3; gm=4/3; nL=100; nB=150
LT = np.linspace(-1,4,nL) #Lambda_T
b0 = np.logspace(-2,3,nB) #beta_0
#A = np.outer(np.ones(nL), g+1 + gm*(gm-g)/(gm+g*b0))
#B = q*( np.outer((g-1)*LT+4-g, np.ones(nB)) 
#+ gm*np.outer(np.ones(nL), 1/b0 - gm*(2-g/gm+1/b0)/(gm+g*b0)) )
#C = q*( np.outer(q*(LT-2)+1, np.ones(nB)) 
#+ q*gm*( np.outer(LT-1,1/b0) + gm*np.outer(np.ones(nL),(1+1/b0)/(gm+g*b0)) ) )

A = np.outer( np.ones(nL), g+1 + gm*(gm+1)/(g*b0) )

B = q*( np.outer((g-1)*LT+4-g, np.ones(nB)) 
- (gm/g)*(np.outer(np.ones(nL), (2*gm-g-4)/b0) - np.outer(LT,(g-1)/b0)) )

C = q*( np.outer(q*(LT-2)+1, np.ones(nB)) 
- (gm/g)*( np.outer( q*(2+g-LT*(1+g)-gm)-1,1/b0) - np.outer(LT,q*gm/b0**2)) )

Det = B*B - 4*A*C
cs = plt.contourf(b0,LT,Det,40,vmin=-1e8,vmax=.1e8)
plt.text(0.015, 3, r'no transonic solution', fontsize=15)
plt.colorbar()
plt.contour(cs,colors='white',levels=[0])
#plt.contour(b0,LT,Det)
plt.xlabel(r'$\beta_0$')
plt.ylabel(r'$\Lambda_T$')
plt.xscale('log')
plt.title(r'$\Delta=B^2-4AC$ for $q=2,~\gamma_m=4/3$',fontsize=16)
plt.show()