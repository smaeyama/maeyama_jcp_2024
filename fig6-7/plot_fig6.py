#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

plt.style.use('../nature_style.txt')

ds_phiz = xr.open_dataset("data/rotating_phi_zt.nc")

s_hat = 0.8
gamma_e = 0.2
T_lap = 2*np.pi*s_hat/gamma_e
print(T_lap)
ds_phiz = ds_phiz.interp(t=[T_lap,1.125*T_lap,1.25*T_lap,1.5*T_lap,1.75*T_lap,2*T_lap])
tt = np.array(ds_phiz.t)
zz = np.array(ds_phiz.zz)
rephi = np.array(ds_phiz.rephi)
imphi = np.array(ds_phiz.imphi)
phi = rephi + 1j * imphi
print(ds_phiz)

print(tt[-1]-tt[0],T_lap)
phi0phi=np.sum(np.conjugate(phi[0,:])*phi[-1,:])
phi0_norm=np.sum(np.conjugate(phi[0,:])*phi[0,:])
omega = np.log(phi0phi/phi0_norm)/(-1j*T_lap)
print("Frequency:",omega.real,", growthrate:",omega.imag)
print("Note: Real frequency should have arbitrariness by 2*n*pi/T_lap={:}n.".format(2*np.pi/T_lap))

# labels = ["(a)","(b)","(c)","(d)","(e)","(f)"]
# labels = [r"(a) $\bar{t}=0$",r"(b) $\bar{t}=T_\mathrm{lap}/8$",r"(c) $\bar{t}=T_\mathrm{lap}/4$",r"(d) $\bar{t}=T_\mathrm{lap}/2$",r"(e) $\bar{t}=3T_\mathrm{lap}/4$",r"(d) $\bar{t}=T_\mathrm{lap}$"]
labels = [r"(a) $t=T_\mathrm{lap}$",r"(b) $t=9T_\mathrm{lap}/8$",r"(c) $t=5T_\mathrm{lap}/4$",r"(d) $t=3T_\mathrm{lap}/2$",r"(e) $t=7T_\mathrm{lap}/4$",r"(d) $t=2T_\mathrm{lap}$"] 

# xticks = np.linspace(-3*np.pi,3*np.pi,7,endpoint=True)
# xticklabels = [r"-3$\pi$",r"-2$\pi$",r"-$\pi$","0",r"$\pi$",r"2$\pi$",r"3$\pi$"]
xticks = np.linspace(-2*np.pi,2*np.pi,5,endpoint=True)
xticklabels = [r"-2$\pi$",r"-$\pi$","0",r"$\pi$",r"2$\pi$"]

#Normalization
phi = phi/phi[0,len(zz)//2]


fig=plt.figure(figsize=(7,3.5),dpi=600) # figsize=(width,height(inch)),dpi(dots per inch)
it=0
ax=fig.add_subplot(2,3,it+1)
# ax.set_xlabel(r"Field-aligned coord. $z$")
ax.set_ylabel(r"$\phi_k$ [$T_\mathrm{i}\rho_\mathrm{ti}/(eR_\mathrm{a})$]")
ax.plot(zz,phi[it,:].real,"-",label=r"Re[$\phi$]")
ax.plot(zz,phi[it,:].imag,"--",label=r"Im[$\phi$]")
# ax.set_title(labels[it]+" t={:5.2f}".format(tt[it]))
ax.set_title(labels[it])
ax.set_xlim(-2.5*np.pi,2.5*np.pi)
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
# ax.legend()
for it in range(1,3):
    ax=fig.add_subplot(2,3,it+1)
    # ax.set_xlabel(r"Field-aligned coord. $z$")
    # ax.set_ylabel(r"Poloidal $y$ [$\rho_\mathrm{ti}$]")
    ax.plot(zz,phi[it,:].real,"-",label=r"Re[$\phi_k(z,t)$]")
    ax.plot(zz,phi[it,:].imag,"--",label=r"Im[$\phi_k(z,t)$]")
    # ax.set_title(labels[it]+" t={:5.2f}".format(tt[it]))
    ax.set_title(labels[it])
    ax.set_xlim(-2.5*np.pi,2.5*np.pi)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    # if it==2:
    #     ax.legend(loc="center left",bbox_to_anchor=(1.05,0.5))
it=3
ax=fig.add_subplot(2,3,it+1)
ax.set_xlabel(r"Field-aligned coord. $z''$")
ax.set_ylabel(r"$\phi_k$ [$T_\mathrm{i}\rho_\mathrm{ti}/(eR_\mathrm{a})$]")
ax.plot(zz,phi[it,:].real,"-",label=r"Re[$\phi$]")
ax.plot(zz,phi[it,:].imag,"--",label=r"Im[$\phi$]")
# ax.set_title(labels[it]+" t={:5.2f}".format(tt[it]))
ax.set_title(labels[it])
ax.set_xlim(-2.5*np.pi,2.5*np.pi)
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
# ax.legend()
for it in range(4,6):
    ax=fig.add_subplot(2,3,it+1)
    ax.set_xlabel(r"Field-aligned coord. $z''$")
    # ax.set_ylabel(r"Poloidal $y$ [$\rho_\mathrm{ti}$]")
    ax.plot(zz,phi[it,:].real,"-",label=r"Re[$\phi(t'',z'')$]")
    ax.plot(zz,phi[it,:].imag,"--",label=r"Im[$\phi(t'',z'')$]")
    if it==5:
        ax.plot(zz,(phi[0,:]*np.exp(-1j*omega*T_lap)).real,"+",c="C02",label=r"Re[$\phi_k(t''-T_\mathrm{lap},z'')e^{-i\omega T_\mathrm{lap}}$]",markersize=4)
        ax.plot(zz,(phi[0,:]*np.exp(-1j*omega*T_lap)).imag,"x",c="C03",label=r"Im[$\phi_k(t''-T_\mathrm{lap},z'')e^{-i\omega T_\mathrm{lap}}$]",markersize=3)
    # ax.set_title(labels[it]+" t={:5.2f}".format(tt[it]))
    ax.set_title(labels[it])
    ax.set_xlim(-2.5*np.pi,2.5*np.pi)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    if it==5:
        ax.legend(loc="center left",bbox_to_anchor=(1.05,0.5))

fig.tight_layout()
plt.savefig("fig6.pdf",dpi=600,bbox_inches="tight")
plt.show()


# In[ ]:




