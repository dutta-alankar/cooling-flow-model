#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 12:06:39 2019

@author: alankar

pressure-driven cooling flow in cylindrical, spherical geometry

12 Jan 2021: added magnetic fields
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import solve_ivp
import sys
from scipy import interpolate
import sys
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

def deriv(x, y):
    if y[0]>0:
        print ("negative density!",x)
        sys.exit()
    d = -1./(y[0]*x**q); p = y[1]*d**gamma; # tilde density and pressure 
    T = mu*mp*c0*c0*p/(kB*d*gamma) # temperature in CGS
    Lam = Lambda(T)/Lam0 #tilde Lambda
    Num1 = q*(c0/v0)**2*( d*Lam*(1+gamma_m/(gamma*beta0))/y[0] + (p+gamma_m*d**gamma_m/(gamma*beta0))/(x*d) )
    Den1 =  (1 - (c0/v0)**2*(p+gamma_m*d**gamma_m/(gamma*beta0))/(d*y[0]**2))*y[0]
    return [ Num1/Den1, -q*gamma*(1+gamma_m/(gamma*beta0))*Lam*d**(2-gamma)/y[0] ]

def Lambda(temp): #returns cooling function in cgs
    klo=0; khi=tab_sz-1
    while (klo != (khi-1)):
        kmid = int((khi+klo)/2)
        Tmid = Ttab[kmid]
        if (temp<=Tmid):
            khi = kmid
        if (temp>Tmid):
            klo = kmid
    dT = Ttab[khi] - Ttab[klo]
    scrh = Ltab[klo]*(Ttab[khi]-temp)/dT + Ltab[khi]*(temp-Ttab[klo])/dT; #linear interpolation
    return scrh
'''
def Lambda(temp):
    lam = 0.0  
    if (temp<=1.e7 and temp>=1.e4):
        lam = 2.7e-23*(temp/1.e7)**-1.0
    if (temp<1.e4):
        lam = 2.7e-20*(temp/1.e4)**20
        #lam = 1.e-100
    if (temp>1.e7):
        lam = 2.7e-23*(temp/1.e7)**-20
        #lam = 1.e-100
    return lam
'''

D = np.loadtxt('./cooltable.dat')
global Ttab, Ltab, tab_sz
Ttab = D[:,0]; Ltab = D[:,1]; tab_sz = np.size(Ttab)

global q, gamma, gamma_m, beta0, c0, v0, Lam0
gamma_m = 1.03; beta0 = 0.03
q=2; gamma=5./3.
mu = 0.62; mue = 1.17; mui = 1./(1./mu - 1./mue)
v0byc0 = 0.2; T0 = 2.e3; d0 = 1e-2*mu*mp ; c0 = np.sqrt(gamma*kB*T0/mu/mp)
shift = 3.0

ct0=c0*np.sqrt(1+gamma_m/(gamma*beta0)) #sound speed including B-fields
v0 = v0byc0*c0; Lam0 = Lambda(T0); v0byct0 = v0/ct0
tcool0cgs = mp*mue*mui*kB*T0/((gamma-1)*mu*d0*Lam0)
R0 = q*gamma*v0*tcool0cgs*(1+gamma_m/(gamma*beta0))

if (v0byct0==1):
    epsT=0.01
    LamT = np.log(Lambda(T0*(1+epsT))/Lambda(T0*(1-epsT)))/np.log((1+epsT)/(1-epsT)) #dlnLambda/dlnT
    Aq = (gamma+1) + gamma_m*(1+gamma_m)/(gamma*beta0)
    Bq = q*( LamT*(gamma-1)+4-gamma -gamma_m*(2*gamma_m-gamma-4-LamT*(gamma-1))/(gamma*beta0) )
    Cq = q*( q*(LamT-2)+1 - gamma_m*(q*(2+gamma-LamT*(1+gamma+gamma_m/beta0)-gamma_m)-1)/(gamma*beta0) )
    Disc = Bq*Bq - 4*Aq*Cq #Discriminant 
    if (Disc<0):
        print ("no transonic solution exists")
        sys.exit()
    vp = (-Bq + np.sqrt(Disc))/(2*Aq) # physically relevant root for the slope of velocity at sonic point
    #eps = 0.00001; 
else:
    #eps = 0.0
    vp = 0.0
eps = 0.00001
sp = (gamma+gamma_m/beta0)*q #slope of entropy at sonic point, same at all critical point
method_ODE = "DOP853"
#method_ODE = "LSODA"
reltol = 5e-13; abstol = 1.e-12
#start at the sonic point and move out
Rout = 5/(R0/kpc); nR = 1000
R_out = np.logspace(np.log10(1.+eps),np.log10(Rout),nR)
sol_out = solve_ivp(deriv, [1.+eps,R_out[-1]], [-1.+vp*eps,1.+sp*eps], t_eval=R_out, 
                    method=method_ODE, rtol=reltol, atol=abstol) #outer solution
v_out = sol_out.y[0]; d_out = -1./(sol_out.y[0]*sol_out.t**q); p_out = sol_out.y[1]*d_out**gamma
#inner solution
Rin = 0.57 #starting of the inner radius
R = np.logspace(np.log10(1.-eps),np.log10(Rin),nR)
sol = solve_ivp(deriv, [1.-eps,R[-1]], [-1.-vp*eps, 1.-sp*eps], t_eval=R, 
                method=method_ODE, rtol=reltol, atol=abstol) #inner solution
v = sol.y[0]; d = -1./(sol.y[0]*sol.t**q); p = sol.y[1]*d**gamma

#analyze results
R = np.concatenate((np.flip(R), R_out)); v = np.concatenate((np.flip(v), v_out)) #arranging inner & outer solutions in increasing order
d = np.concatenate((np.flip(d), d_out)); p = np.concatenate((np.flip(p), p_out))
Mach = -v*(v0/c0)/np.sqrt( p/d + gamma_m*d**(gamma_m-1)/(beta0*gamma) ) #including B-fields in sound speed
Be = 0.5*v**2 + p*(c0/v0)**2/((gamma-1.)*d) + gamma_m*d**(gamma_m-1)*(c0/v0)**2/(gamma*(gamma_m-1)*beta0) #including B-fields, normalized to v_0^2

T = mu*mp*c0*c0*p/(gamma*kB*d); Lam=np.zeros(np.size(T))
for i in range(np.size(T)):
    Lam[i]=Lambda(T[i])/Lam0

vcgs = v0*v 
pcgs = p*d0*c0**2/gamma
pmagcgs = d0*c0**2*d**gamma_m/(gamma*beta0)
beta = pcgs/pmagcgs
Rcgs = R0*R
dcgs = d*d0
vdvdr = vcgs*np.gradient(vcgs,Rcgs)
dpdrbyrho = np.gradient(pcgs,Rcgs)/dcgs
dpmagbydrho = np.gradient(pmagcgs,Rcgs)/dcgs
test_mom = vdvdr + dpdrbyrho + dpmagbydrho
test_mom /= vdvdr
dsbydr = np.gradient(np.log(pcgs/dcgs**gamma),Rcgs)
Tcgs = (mu*mp/kB)*pcgs/dcgs
pmag_cgs = (d0*c0**2/gamma)*(dcgs/d0)**gamma_m/beta0
neniL = dcgs**2*Lam*Lam0/(mue*mui*mp**2)
test_energy = pcgs*vcgs*dsbydr/(gamma-1) + neniL 
test_energy = test_energy/neniL
Becgs = Be*v0**2
Mdot = -Rcgs**q*dcgs*vcgs

checks = False
if checks:
    if (q==2):
        Mdot *= 4.*np.pi
        dEdotdT = 4.*np.pi*Rcgs**2*neniL/np.gradient(Tcgs,Rcgs)
        MdotdBedT = Mdot*np.gradient(Becgs,Tcgs)
        plt.loglog(Tcgs,Tcgs*dEdotdT,label="dEdot/dT")
        plt.loglog(Tcgs,Tcgs*MdotdBedT,label="Mdot*dBe/dT")
        plt.ylabel(r'erg s$^{-1}$ K$^{-1}$')
        plt.xlabel('T(K)')
        plt.legend()
    plt.figure()
    plt.plot(R,test_mom,label="momentum error; [vdv/dr+dp/dr/rho+dpmag/dr/rho]/vdv/dr")
    plt.plot(R,test_energy,label="energy error; [pvds/dr/(gamma-1) + neniL]/neniL")
    plt.legend()
    plt.xlabel('dedimensionalized radius')
    plt.figure()
    
    plt.loglog(R,-v,label='velocity')
    plt.loglog(R,p,label='pressure')
    plt.loglog(R,d,label='density')
    plt.loglog(R,T,label='temperature')
    plt.loglog(R,p/d**gamma,label='entropy')
    plt.xlabel('dedimensionalized radius')
    plt.grid()
    plt.legend()
    plt.show()

np.savetxt('bernoulli-model(rTBeMdot).txt',
           np.vstack((Rcgs,Tcgs,Becgs,Mdot)).T)

#--------------------- Comparing ---------------------------

sizes = [0.5, 1.0, 1.5, 2.0] #kpc
select = 2
r_crit = sizes[select]
start = select*3+1

#----------------------------------
dn_data = np.loadtxt('./dylan-nelson/fig10_TNG50-1_z0.5_h8_nH.txt', skiprows=5)
dn_radius = dn_data[:,0]
#print(dn_radius)
dn_nH_50 = dn_data[:,start]
dn_nH_16, dn_nH_84 = dn_data[:,start+1],dn_data[:,start+2]
profile_s4 = interpolate.interp1d(dn_radius, dn_nH_50, fill_value='extrapolate')
nH_crit = 10**profile_s4(r_crit) #cm^-3

#----------------------------------
dn_data = np.loadtxt('./dylan-nelson/fig10_TNG50-1_z0.5_h8_beta.txt', skiprows=5)
dn_radius = dn_data[:,0]
#print(dn_radius)
dn_beta_50 = dn_data[:,start]
dn_beta_16, dn_beta_84 = dn_data[:,start+1],dn_data[:,start+2]
profile_s4 = interpolate.interp1d(dn_radius, dn_beta_50, fill_value='extrapolate')
beta_crit = 10**profile_s4(r_crit) #dimensionless
#beta_crit = 10**np.min(dn_beta_s4)

#----------------------------------
dn_data = np.loadtxt('./dylan-nelson/fig10_TNG50-1_z0.5_h8_temp.txt', skiprows=5)
dn_radius = dn_data[:,0]
#print(dn_radius)
dn_temp_50 = dn_data[:,start]
dn_temp_16, dn_temp_84 = dn_data[:,start+1],dn_data[:,start+2]
profile_s4 = interpolate.interp1d(dn_radius, dn_temp_16, fill_value='extrapolate')
temp_crit = 10**profile_s4(r_crit) #K

#----------------------------------
dn_data = np.loadtxt('./dylan-nelson/fig10_TNG50-1_z0.5_h8_vel_rel.txt', skiprows=5)
dn_radius = dn_data[:,0]
#print(dn_radius)
dn_vel_50 = dn_data[:,start]
dn_vel_16, dn_vel_84 = dn_data[:,start+1],dn_data[:,start+2]
profile_s4 = interpolate.interp1d(dn_radius, dn_vel_50, fill_value='extrapolate')
vel_crit = profile_s4(r_crit)*1e5 #cm s^-1
#vel_crit = np.min(dn_vel_s4)*1e5

#----------------------------------
dn_data = np.loadtxt('./dylan-nelson/fig10_TNG50-1_z0.5_h8_P_tot.txt', skiprows=5)
dn_radius = dn_data[:,0]
#print(dn_radius)
dn_Ptot_50 = dn_data[:,start]
dn_Ptot_16, dn_Ptot_84 = dn_data[:,start+1],dn_data[:,start+2]
profile_s4 = interpolate.interp1d(dn_radius, dn_Ptot_50, fill_value='extrapolate')
Ptot_crit = 10**profile_s4(r_crit) #K cm^-3

#----------------------------------
dn_data = np.loadtxt('./dylan-nelson/fig10_TNG50-1_z0.5_h8_P_gas.txt', skiprows=5)
dn_radius = dn_data[:,0]
#print(dn_radius)
dn_Pgas_50 = dn_data[:,start]
dn_Pgas_16, dn_Pgas_84 = dn_data[:,start+1],dn_data[:,start+2]
profile_s4 = interpolate.interp1d(dn_radius, dn_Pgas_50, fill_value='extrapolate')
Pgas_crit = 10**profile_s4(r_crit) #K cm^-3

#----------------------------------

Rcgs *= shift

#Plotting model and Dylan data together
fig, axs = plt.subplots(2, 3, figsize=(20,11))
gs1 = gridspec.GridSpec(2, 3)
gs1.update(wspace=0.28, hspace=0.02) # set the spacing between axes.
axs = np.array([ [plt.subplot(gs1[0]), plt.subplot(gs1[1]), plt.subplot(gs1[2])], 
                 [plt.subplot(gs1[3]), plt.subplot(gs1[4]), plt.subplot(gs1[5])] ])

axs[0, 0].plot(Rcgs/kpc, np.log10(dcgs/(mu*mp)), linewidth=5)
axs[0, 0].plot(dn_radius, dn_nH_50, c='tab:gray', linewidth=5)
axs[0, 0].fill_between(dn_radius, dn_nH_16, dn_nH_84, color='tab:gray', alpha=0.3)
axs[0, 0].grid()
axs[0, 0].set_ylabel(r'number density [$\rm cm^{-3}\ (log)$]', size=20 )
axs[0, 0].get_xaxis().set_ticklabels([])
axs[0, 0].set_xlim(xmin=0., xmax=5.0)
axs[0, 0].tick_params(axis='both', which='major', labelsize=18)
axs[0, 0].tick_params(axis='both', which='minor', labelsize=16)
axs[0, 0].xaxis.set_major_locator(plt.MaxNLocator(6))

axs[0, 1].plot(Rcgs/kpc, np.log10(pcgs/kB), 'tab:orange', linewidth=5)
axs[0, 1].plot(dn_radius, dn_Pgas_50, c='tab:gray', linewidth=5)
axs[0, 1].fill_between(dn_radius, dn_Pgas_16, dn_Pgas_84, color='tab:gray', alpha=0.3)
axs[0, 1].grid()
axs[0, 1].set_ylim(ymin=1.8)
axs[0, 1].set_xlim(xmin=0., xmax=5.0)
axs[0, 1].set_ylabel(r'$\rm p_{gas}$ [$\rm \times k_B \ K cm^{-3}\ (log)$]', size=20 )
axs[0, 1].get_xaxis().set_ticklabels([])
axs[0, 1].tick_params(axis='both', which='major', labelsize=18)
axs[0, 1].tick_params(axis='both', which='minor', labelsize=16)
axs[0, 1].xaxis.set_major_locator(plt.MaxNLocator(6))

axs[0, 2].plot(Rcgs/kpc, vcgs/1e5, 'tab:green', linewidth=5)
axs[0, 2].plot(dn_radius, dn_vel_50, c='tab:gray', linewidth=5)
axs[0, 2].fill_between(dn_radius, dn_vel_16, dn_vel_84, color='tab:gray', alpha=0.3)
axs[0, 2].grid()
axs[0, 2].set_ylabel(r'velocity [$\rm km s^{-1}$]', size=20 )
axs[0, 2].set_xlim(xmin=0., xmax=5.0)
axs[0, 2].get_xaxis().set_ticklabels([])
axs[0, 2].tick_params(axis='both', which='major', labelsize=18)
axs[0, 2].tick_params(axis='both', which='minor', labelsize=16)
axs[0, 2].xaxis.set_major_locator(plt.MaxNLocator(6))

axs[1, 0].plot(Rcgs/kpc, np.log10(beta), 'tab:red', linewidth=5)
axs[1, 0].plot(dn_radius, dn_beta_50, c='tab:gray', linewidth=5)
axs[1, 0].fill_between(dn_radius, dn_beta_16, dn_beta_84, color='tab:gray', alpha=0.3)
axs[1, 0].grid()
axs[1, 0].set_xlabel(r'distance [$\rm kpc$]', size=20 )
axs[1, 0].set_ylabel(r'$\rm \beta = p_{gas}/p_{mag}\ (log)$', size=20 )
axs[1, 0].set_xlim(xmin=0., xmax=5.0)
axs[1, 0].set_ylim(ymin=-1.2, ymax=2.2)
axs[1, 0].tick_params(axis='both', which='major', labelsize=18)
axs[1, 0].tick_params(axis='both', which='minor', labelsize=16)
axs[1, 0].xaxis.set_major_locator(plt.MaxNLocator(6))

axs[1, 1].plot(Rcgs/kpc, np.log10(Tcgs), 'tab:brown', linewidth=5)
axs[1, 1].plot(dn_radius, dn_temp_50, c='tab:gray', linewidth=5)
axs[1, 1].fill_between(dn_radius, dn_temp_16, dn_temp_84, color='tab:gray', alpha=0.3)
axs[1, 1].grid()
axs[1, 1].set_xlim(xmin=0., xmax=5.0)
axs[1, 1].set_ylim(ymin=3.5)
axs[1, 1].set_ylabel(r' temperature [$\rm K\ (log)$]', size=20 )
axs[1, 1].set_xlabel(r'distance [$\rm kpc$]', size=20 )
axs[1, 1].tick_params(axis='both', which='major', labelsize=18)
axs[1, 1].tick_params(axis='both', which='minor', labelsize=16)
axs[1, 1].xaxis.set_major_locator(plt.MaxNLocator(6))

axs[1, 2].plot(Rcgs/kpc, np.log10((pcgs+pmagcgs)/kB), 'tab:purple', linewidth=5)
axs[1, 2].plot(dn_radius, dn_Ptot_50, c='tab:gray', linewidth=5)
axs[1, 2].fill_between(dn_radius, dn_Ptot_16, dn_Ptot_84, color='tab:gray', alpha=0.3)
axs[1, 2].grid()
axs[1, 2].set_xlim(xmin=0., xmax=5.0)
axs[1, 2].set_ylabel(r'$\rm p_{total}$ [$\rm \times k_B \ K cm^{-3}\ (log)$] ', size=20 )
axs[1, 2].set_xlabel(r'distance [$\rm kpc$]', size=20 )
axs[1, 2].tick_params(axis='both', which='major', labelsize=18)
axs[1, 2].tick_params(axis='both', which='minor', labelsize=16)
axs[1, 2].xaxis.set_major_locator(plt.MaxNLocator(6))
'''
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
'''
Mdot_avg = np.mean(Mdot/(Msun/yr))
#fig.suptitle('Critical radius is approximately set to %.2f kpc\n'%(R0/kpc*shift)+r'$\dot{\rm M} \rm \approx %.1f \times 10^{%d}\ M_{\odot} \ yr^{-1}$'%( fman(Mdot_avg) ,fexp(Mdot_avg)), size=25)
fig.subplots_adjust(top=0.92)
plt.savefig('example-prof-marchinout.png', transparent=True, bbox_inches='tight')
plt.show()
