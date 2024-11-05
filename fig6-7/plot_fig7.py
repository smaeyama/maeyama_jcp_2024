#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

# plt.style.use('../nature_style.txt')

ds_phiz = xr.open_dataset("data/rotating_phi_zt.nc")
tt = np.array(ds_phiz.t)
zz = np.array(ds_phiz.zz)
rephi = np.array(ds_phiz.rephi)
imphi = np.array(ds_phiz.imphi)
phi = rephi + 1j * imphi


s_hat = 0.8
gamma_e = 0.2
T_lap = 2*np.pi*s_hat/gamma_e
temp_phi = ds_phiz.interp(t=[T_lap,2*T_lap])
temp_phi = np.array(temp_phi.rephi) + 1j * np.array(temp_phi.imphi)
phi0phi=np.sum(np.conjugate(temp_phi[0,:])*temp_phi[-1,:])
phi0_norm=np.sum(np.conjugate(temp_phi[0,:])*temp_phi[0,:])
omega = np.log(phi0phi/phi0_norm)/(-1j*T_lap)
print("Frequency:",omega.real,", growthrate:",omega.imag)
print("Note: Real frequency should have arbitrariness by 2*n*pi/T_lap={:}n.".format(2*np.pi/T_lap))

amplitude = np.log(phi[:,len(zz)//2]).real
plt.plot(tt,amplitude)
plt.plot(tt,np.log(np.abs(phi[:,len(zz)//2])),".")
plt.show()

phase = np.log(phi[:,len(zz)//2]).imag
plt.plot(tt,phase)
plt.show()

wphase = phase.copy()
for i in range(1,len(phase)):
    if np.abs(wphase[i]-wphase[i-1]) > np.pi:
        if wphase[i]>wphase[i-1]:
            phase[i:] = phase[i:]-2*np.pi
        else:
            phase[i:] = phase[i:]+2*np.pi
plt.plot(tt,phase)
plt.plot(tt,-omega.real*(tt-tt[0])+phase[0],label=r"$\omega_r$={:}".format(omega.real))
plt.plot(tt,-(omega.real-2*np.pi/T_lap)*(tt-tt[0])+phase[0],label=r"$\omega_r$={:}".format(omega.real-2*np.pi/T_lap))
plt.plot(tt,-(omega.real-2*(2*np.pi/T_lap))*(tt-tt[0])+phase[0],label=r"$\omega_r$={:}".format(omega.real-2*(2*np.pi/T_lap)))
plt.plot(tt,-(omega.real-3*(2*np.pi/T_lap))*(tt-tt[0])+phase[0],label=r"$\omega_r$={:}".format(omega.real-3*(2*np.pi/T_lap)))
plt.legend()
plt.show()


# In[ ]:


ds_phiz_gamma0 = xr.open_dataset("data/gamma0_phi_zt.nc")
tt_gamma0 = np.array(ds_phiz_gamma0.t)
zz_gamma0 = np.array(ds_phiz_gamma0.zz)
rephi_gamma0 = np.array(ds_phiz_gamma0.rephi)
imphi_gamma0 = np.array(ds_phiz_gamma0.imphi)
phi_gamma0 = rephi_gamma0 + 1j * imphi_gamma0

omega_gamma0 = np.log(phi_gamma0[-1,len(zz_gamma0)//2]/phi_gamma0[-2,len(zz_gamma0)//2])/(-1j*(tt_gamma0[-1]-tt_gamma0[-2]))
print(r"Frequency($\gamma_E=0$):",omega_gamma0.real,r", growthrate($\gamma_E=0$):",omega_gamma0.imag)


amplitude_gamma0 = np.log(phi_gamma0[:,len(zz_gamma0)//2]).real
plt.plot(tt_gamma0,amplitude_gamma0)
plt.plot(tt_gamma0,np.log(np.abs(phi_gamma0[:,len(zz_gamma0)//2])),".")
plt.plot(tt_gamma0,amplitude_gamma0[0]+omega_gamma0.imag*(tt_gamma0-tt_gamma0[0]))
plt.show()

phase_gamma0 = np.log(phi_gamma0[:,len(zz_gamma0)//2]).imag
plt.plot(tt_gamma0,phase_gamma0)
plt.show()

wphase = phase_gamma0.copy()
for i in range(1,len(phase_gamma0)):
    if np.abs(wphase[i]-wphase[i-1]) > np.pi:
        if wphase[i]>wphase[i-1]:
            phase_gamma0[i:] = phase_gamma0[i:]-2*np.pi
        else:
            phase_gamma0[i:] = phase_gamma0[i:]+2*np.pi
plt.plot(tt_gamma0,phase_gamma0)
plt.plot(tt_gamma0,-omega_gamma0.real*(tt_gamma0-tt_gamma0[0])+phase_gamma0[0],label=r"$\omega_r$={:}".format(omega_gamma0.real))
plt.legend()
plt.show()


# In[ ]:


# import numpy as np
# import matplotlib.pyplot as plt
# import xarray as xr

# plt.style.use('../nature_style.txt')

# ds_phiz = xr.open_dataset("data/rotating_phi_zt.nc")
# tt = np.array(ds_phiz.t)
# zz = np.array(ds_phiz.zz)
# rephi = np.array(ds_phiz.rephi)
# imphi = np.array(ds_phiz.imphi)
# phi = rephi + 1j * imphi


# s_hat = 0.8
# gamma_e = 0.2
# T_lap = 2*np.pi*s_hat/gamma_e
# temp_phi = ds_phiz.interp(t=[T_lap,2*T_lap])
# temp_phi = np.array(temp_phi.rephi) + 1j * np.array(temp_phi.imphi)
# phi0phi=np.sum(np.conjugate(temp_phi[0,:])*temp_phi[-1,:])
# phi0_norm=np.sum(np.conjugate(temp_phi[0,:])*temp_phi[0,:])
# omega = np.log(phi0phi/phi0_norm)/(-1j*T_lap)
# print("Frequency:",omega.real,", growthrate:",omega.imag)
# print("Note: Real frequency should have arbitrariness by 2*n*pi/T_lap={:}n.".format(2*np.pi/T_lap))

# amplitude_theta0 = []
# phase_theta0 = []
# for wt in tt:
#     theta = 0
#     wz = theta - gamma_e * (wt-tt[0]) / s_hat
#     wds = ds_phiz.interp(t=wt,zz=wz)
#     # print(wds)
#     wphi = np.array(wds.rephi + 1j * wds.imphi)
#     wamplitude = np.log(wphi).real
#     wphase = np.log(wphi).imag
#     amplitude_theta0.append(wamplitude)
#     phase_theta0.append(wphase)
# amplitude_theta0 = np.array(amplitude_theta0)
# phase_theta0 = np.array(phase_theta0)
# plt.plot(tt,amplitude)
# plt.plot(tt,amplitude_theta0)
# plt.show()
# plt.plot(tt,phase_theta0)
# plt.show()

# wphase = phase_theta0.copy()
# for i in range(1,len(phase_theta0)):
#     if np.abs(wphase[i]-wphase[i-1]) > np.pi:
#         if wphase[i]>wphase[i-1]:
#             phase_theta0[i:] = phase_theta0[i:]-2*np.pi
#         else:
#             phase_theta0[i:] = phase_theta0[i:]+2*np.pi


# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

plt.style.use('../nature_style.txt')

ds_phiz = xr.open_dataset("data/rotating_phi_zt.nc")
tt = np.array(ds_phiz.t)
zz = np.array(ds_phiz.zz)
rephi = np.array(ds_phiz.rephi)
imphi = np.array(ds_phiz.imphi)
phi = rephi + 1j * imphi

s_hat = 0.8
gamma_e = 0.2
T_lap = 2*np.pi*s_hat/gamma_e
temp_phi = ds_phiz.interp(t=[T_lap,2*T_lap])
temp_phi = np.array(temp_phi.rephi) + 1j * np.array(temp_phi.imphi)
phi0phi=np.sum(np.conjugate(temp_phi[0,:])*temp_phi[-1,:])
phi0_norm=np.sum(np.conjugate(temp_phi[0,:])*temp_phi[0,:])
omega = np.log(phi0phi/phi0_norm)/(-1j*T_lap)
omega = omega-3*(2*np.pi/T_lap)
print("Frequency:",omega.real,", growthrate:",omega.imag)
print("Note: Real frequency should have arbitrariness by 2*n*pi/T_lap={:}n.".format(2*np.pi/T_lap))

tt = np.linspace(T_lap,2*T_lap,int(T_lap/0.1)+1)
ds_phiz_interp = ds_phiz.interp(t=tt)
zz = np.array(ds_phiz_interp.zz)
rephi = np.array(ds_phiz_interp.rephi)
imphi = np.array(ds_phiz_interp.imphi)
phi = rephi + 1j * imphi
phi = phi/phi[0,len(zz)//2]

amplitude = np.log(phi[:,len(zz)//2]).real
phase = np.log(phi[:,len(zz)//2]).imag
wphase = phase.copy()
for i in range(1,len(phase)):
    if np.abs(wphase[i]-wphase[i-1]) > np.pi:
        if wphase[i]>wphase[i-1]:
            phase[i:] = phase[i:]-2*np.pi
        else:
            phase[i:] = phase[i:]+2*np.pi

plt.rcParams["xtick.top"] = False

fig=plt.figure(figsize=(3,4),dpi=600) # figsize=(width,height(inch)),dpi(dots per inch)
ax=fig.add_subplot(2,1,1)
ax.set_xlabel(r"Time $t$ [$R_\mathrm{a}/v_\mathrm{ti}$]")
ax.set_ylabel(r"Amplitude $|\phi_k|(t'',z''=0)$ [a.u.]")
ax.plot(tt,np.exp(amplitude),"-",label=r"|$\phi_k$|")
ax.plot(tt[:],np.exp(amplitude[0]+omega.imag*(tt[:]-tt[0])),"--",lw=1,label=r"$e^{\gamma_l t}$")
ax.plot(tt[:len(tt)//3],np.exp(amplitude[0]+omega_gamma0.imag*(tt[:len(tt)//3]-tt[0])),":",c="k",lw=1,label=r"$e^{\gamma_l^{(\gamma_E=0)} t}$")
ax.set_yscale("log")
ax.set_xlim(T_lap,2.05*T_lap)
ax.set_ylim(1,None)
ax.set_xticks(np.linspace(T_lap,2*T_lap,5))
ax.set_xticklabels([r"$T_\mathrm{lap}$",r"$5T_\mathrm{lap}/4$",r"$3T_\mathrm{lap}/2$",r"$7T_\mathrm{lap}/4$",r"$2T_\mathrm{lap}$"])
t0=tt[0]
def forward(t):
    zz=0
    theta = zz + gamma_e * (t-t0) / s_hat
    return theta
def inverse(theta):
    zz = 0
    t = t0 + (theta - zz) * s_hat / gamma_e
    return t
secax = ax.secondary_xaxis("top", functions=(forward, inverse))
secax.set_xlabel(r"Poloidal angle $\theta(t'',z''=0)=t''\gamma_E/\hat{s}$")
secax.set_xticks(np.linspace(0,2*np.pi,5))
secax.set_xticklabels(["0",r"$\pi/2$",r"$\pi$",r"$3\pi/2$",r"$2\pi$"])
ax.legend(loc="lower right")

ax=fig.add_subplot(2,1,2)
ax.set_xlabel(r"Time $t$ [$R_\mathrm{a}/v_\mathrm{ti}$]")
ax.set_ylabel(r"Phase $\arg \phi_k(t'',z''=0)$")
ax.plot(tt,phase,"-",label=r"$\arg \phi_k$")
ax.plot(tt,-omega.real*(tt[:]-tt[0]),"--",lw=1,label=r"$-\omega_r t$")
ax.plot(tt[:len(tt)//3],-omega_gamma0.real*(tt[:len(tt)//3]-tt[0]),":",c="k",lw=1,label=r"$-\omega_r^{(\gamma_E=0)} t$")
ax.set_xlim(T_lap,2.05*T_lap)
ax.set_xticks(np.linspace(T_lap,2*T_lap,5))
ax.set_xticklabels([r"$T_\mathrm{lap}$",r"$5T_\mathrm{lap}/4$",r"$3T_\mathrm{lap}/2$",r"$7T_\mathrm{lap}/4$",r"$2T_\mathrm{lap}$"])
ax.set_ylim(0,None)
t0=tt[0]
def forward(t):
    zz=0
    theta = zz + gamma_e * (t-t0) / s_hat
    return theta
def inverse(theta):
    zz = 0
    t = t0 + (theta - zz) * s_hat / gamma_e
    return t
secax = ax.secondary_xaxis("top", functions=(forward, inverse))
secax.set_xlabel(r"Poloidal angle $\theta(t'',z''=0)=t''\gamma_E/\hat{s}$")
secax.set_xticks(np.linspace(0,2*np.pi,5))
secax.set_xticklabels(["0",r"$\pi/2$",r"$\pi$",r"$3\pi/2$",r"$2\pi$"])
ax.legend(loc="lower right")

fig.text(0.18,0.865,"(a)",color="k", fontfamily="sans-serif", fontweight="bold", fontsize=8)
fig.text(0.18,0.375,"(b)",color="k", fontfamily="sans-serif", fontweight="bold", fontsize=8)

fig.tight_layout()
plt.savefig("fig7.pdf",dpi=600,bbox_inches="tight")
plt.show()


# In[ ]:




