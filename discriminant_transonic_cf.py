#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 17:51:44 2020

@author: prateek
discriminant for the criterion of a real root
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rc('lines', linewidth=3, color='r')
plt.rcParams.update({'font.size': 18})

gamma = 5/3;

alpha = np.linspace(-3,6,10000)
#previous expression
#Dis = (1 - (gamma-1)*(alpha+1))**2 + 4*q*(alpha*q-1)*(gamma+1)
#my new correct
q = 2
Dis = q*q*(gamma*(alpha-1)+4-alpha)**2 - 4*(gamma+1)*q*(q*(alpha-2)+1)
#vp1 = (-q*(gamma*(alpha-1)+4-alpha) + np.sqrt(Dis))/(2*(gamma+1))
#vp2 = (-q*(gamma*(alpha-1)+4-alpha) - np.sqrt(Dis))/(2*(gamma+1))
#plt.plot(alpha,vp1,color='#1f77b4',label='q=2, spherical geometry')
#plt.plot(alpha,vp2,color='#1f77b4')
#Alankar's
#Dis = q*q*(alpha*(gamma-1)+2-gamma)**2 - 4*(gamma+1)*q*(q*(2-alpha)+1)
plt.plot(alpha, Dis, '-', label='q=2, spherical geometry')
q = 1
Dis = q*q*(gamma*(alpha-1)+4-alpha)**2 - 4*(gamma+1)*q*(q*(alpha-2)+1)
#vp1 = (-q*(gamma*(alpha-1)+4-alpha) + np.sqrt(Dis))/(2*(gamma+1))
#vp2 = (-q*(gamma*(alpha-1)+4-alpha) - np.sqrt(Dis))/(2*(gamma+1))
#plt.plot(alpha,vp1,'--', color='#ff7f0e',label='q=1, cylindrical geometry')
#plt.plot(alpha,vp2,'--',color='#ff7f0e')
#plt.grid()
#plt.show()
#plt.xlim([-2,5])
#plt.ylim([-3,2])
#plt.legend()
#plt.xlabel(r'$\Lambda_T$')
#plt.ylabel(r'two roots of $\tilde{v}^\prime$ at sonic point', fontsize=16)
#plt.fill([-2, 5, 5, -2],[-3, -3, 0, 0],alpha=0.2)
#plt.text(-1.5, -1.8, r'negative roots', fontsize=16)
plt.plot(alpha, Dis, '--', label='q=1, cylindrical geometry')
plt.legend()
plt.grid()
plt.xlim([-2,5])
plt.ylim([-25,150])
plt.xlabel(r'$\Lambda_T$')
plt.ylabel(r'Discriminant ($\Delta$)')
plt.fill([-2, 5, 5, -2],[-25, -25, 0, 0],alpha=0.2)
plt.text(-1, -14, r'$\Delta$<0, no transonic solution', fontsize=16)