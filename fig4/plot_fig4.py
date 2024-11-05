#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

plt.style.use('../nature_style.txt')

ds_phixy = xr.open_dataset("data/remap_phi_xyt_z0.nc")
ds_phikx = xr.open_dataset("data/remap_phi_kxt_ky05.nc")
tt = np.array(ds_phixy.t)
xx = np.array(ds_phixy.xx)
yy = np.array(ds_phixy.yy)
kx = np.array(ds_phikx.kx)
# labels = ["(a)","(b)","(c)","(d)"]
# labels = [r"(a) $\bar{t}=0$",r"(b) $\bar{t}=T_\mathrm{remap}/2$",r"(c) $\bar{t}=T_\mathrm{remap}$",r"(d) $\bar{t}=3T_\mathrm{remap}/2$"]
labels = [r"(a) $t=T_\mathrm{lap}$",r"(b) $t=T_\mathrm{lap}+T_\mathrm{remap}/2$",r"(c) $t=T_\mathrm{lap}+T_\mathrm{remap}$",r"(d) $t=T_\mathrm{lap}+3T_\mathrm{remap}/2$"]

fig=plt.figure(figsize=(7,2.8),dpi=600) # figsize=(width,height(inch)),dpi(dots per inch)
gs=fig.add_gridspec(nrows=2,ncols=4,height_ratios=(2,1))
it=0
ax=fig.add_subplot(gs[0,0])
ax.set_xlabel(r"Radial $x$ [$\rho_\mathrm{ti}$]")
ax.set_ylabel(r"Poloidal $y$ [$\rho_\mathrm{ti}$]")
ax.pcolormesh(xx,yy,ds_phixy.phi[it,:,:],cmap="jet", rasterized=True)
ax.set_aspect("equal")
# ax.set_title(labels[it]+" t={:5.2f}".format(tt[it]))
ax.set_title(labels[it])
for it in range(1,len(tt)):
    ax=fig.add_subplot(gs[0,it])
    ax.set_xlabel(r"Radial $x$ [$\rho_\mathrm{ti}$]")
    # ax.set_ylabel(r"Poloidal $y$ [$\rho_\mathrm{ti}$]")
    ax.pcolormesh(xx,yy,ds_phixy.phi[it,:,:],cmap="jet", rasterized=True)
    ax.set_aspect("equal")
    # ax.set_title(labels[it]+" t={:5.2f}".format(tt[it]))
    ax.set_title(labels[it])

# labels = ["(e)","(f)","(g)","(h)"]
# labels = [r"(e) $\bar{t}=T_\mathrm{0}$",r"(f) $\bar{t}=T_\mathrm{remap}/2$",r"(g) $\bar{t}=T_\mathrm{remap}$",r"(h) $\bar{t}=3T_\mathrm{remap}/2$"]
labels = [r"(e) $t=T_\mathrm{lap}$",r"(f) $t=T_\mathrm{lap}+T_\mathrm{remap}/2$",r"(g) $t=T_\mathrm{lap}+T_\mathrm{remap}$",r"(h) $t=T_\mathrm{lap}+3T_\mathrm{remap}/2$"]
it=0
ax=fig.add_subplot(gs[1,0])
ax.set_xlabel(r"Wavenumber $k_x$ [$\rho_\mathrm{ti}^{-1}$]")
ax.set_ylabel(r"$\sqrt{\langle|\phi|^2\rangle}$ [$T_\mathrm{i}\rho_\mathrm{ti}/(eR_a)$]")
ax.plot(kx,np.sqrt(ds_phikx.phi[it,:]),".-")
ax.set_xlim(-4,4)
ax.set_ylim(0,None)
# ax.set_title(labels[it]+" t={:5.2f}".format(tt[it]))
ax.set_title(labels[it])
for it in range(1,len(tt)):
    ax=fig.add_subplot(gs[1,it])
    ax.set_xlabel(r"Wavenumber $k_x$ [$\rho_\mathrm{ti}^{-1}$]")
    # ax.set_ylabel(r"$\langle|\phi|^2\rangle$ [$T_\mathrm{i}\rho_\mathrm{ti}/(eR_a)$]")
    ax.plot(kx,np.sqrt(ds_phikx.phi[it,:]),".-")
    ax.set_xlim(-4,4)
    ax.set_ylim(0,None)
    # ax.set_title(labels[it]+" t={:5.2f}".format(tt[it]))
    ax.set_title(labels[it])

fig.tight_layout()
plt.savefig("fig4_remap.pdf",dpi=600,bbox_inches="tight")
plt.show()


# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

plt.style.use('../nature_style.txt')

ds_phixy = xr.open_dataset("data/rotating_phi_xyt_z0.nc")
ds_phikx = xr.open_dataset("data/rotating_phi_kxt_ky05.nc")
tt = np.array(ds_phixy.t)
xx = np.array(ds_phixy.xx)
yy = np.array(ds_phixy.yy)
kx = np.array(ds_phikx.kx)
# labels = ["(i)","(j)","(k)","(l)"]
# labels = [r"(i) $\bar{t}=0$",r"(j) $\bar{t}=T_\mathrm{remap}/2$",r"(k) $\bar{t}=T_\mathrm{remap}$",r"(l) $\bar{t}=3T_\mathrm{remap}/2$"]
labels = [r"(i) $t=T_\mathrm{lap}$",r"(j) $t=T_\mathrm{lap}+T_\mathrm{remap}/2$",r"(k) $t=T_\mathrm{lap}+T_\mathrm{remap}$",r"(l) $t=T_\mathrm{lap}+3T_\mathrm{remap}/2$"]

fig=plt.figure(figsize=(7,2.8),dpi=600) # figsize=(width,height(inch)),dpi(dots per inch)
gs=fig.add_gridspec(nrows=2,ncols=4,height_ratios=(2,1))
it=0
ax=fig.add_subplot(gs[0,0])
ax.set_xlabel(r"Radial $x$ [$\rho_\mathrm{ti}$]")
ax.set_ylabel(r"Poloidal $y$ [$\rho_\mathrm{ti}$]")
ax.pcolormesh(xx,yy,ds_phixy.phi[it,:,:],cmap="jet", rasterized=True)
ax.set_aspect("equal")
# ax.set_title(labels[it]+" t={:5.2f}".format(tt[it]))
ax.set_title(labels[it])
for it in range(1,len(tt)):
    ax=fig.add_subplot(gs[0,it])
    ax.set_xlabel(r"Radial $x$ [$\rho_\mathrm{ti}$]")
    # ax.set_ylabel(r"Poloidal $y$ [$\rho_\mathrm{ti}$]")
    ax.pcolormesh(xx,yy,ds_phixy.phi[it,:,:],cmap="jet", rasterized=True)
    ax.set_aspect("equal")
    # ax.set_title(labels[it]+" t={:5.2f}".format(tt[it]))
    ax.set_title(labels[it])

# labels = ["(m)","(n)","(o)","(p)"]
# labels = [r"(m) $\bar{t}=0$",r"(n) $\bar{t}=T_\mathrm{remap}/2$",r"(o) $\bar{t}=T_\mathrm{remap}$",r"(p) $\bar{t}=3T_\mathrm{remap}/2$"]
labels = [r"(m) $t=T_\mathrm{lap}$",r"(n) $t=T_\mathrm{lap}+T_\mathrm{remap}/2$",r"(o) $t=T_\mathrm{lap}+T_\mathrm{remap}$",r"(p) $t=T_\mathrm{lap}+3T_\mathrm{remap}/2$"]
it=0
ax=fig.add_subplot(gs[1,0])
ax.set_xlabel(r"Wavenumber $k_x$ [$\rho_\mathrm{ti}^{-1}$]")
ax.set_ylabel(r"$\sqrt{\langle|\phi|^2\rangle}$ [$T_\mathrm{i}\rho_\mathrm{ti}/(eR_a)$]")
ax.plot(kx,np.sqrt(ds_phikx.phi[it,:]),".-")
ax.set_xlim(-4,4)
ax.set_ylim(0,None)
# ax.set_title(labels[it]+" t={:5.2f}".format(tt[it]))
ax.set_title(labels[it])
for it in range(1,len(tt)):
    ax=fig.add_subplot(gs[1,it])
    ax.set_xlabel(r"Wavenumber $k_x$ [$\rho_\mathrm{ti}^{-1}$]")
    # ax.set_ylabel(r"$\langle|\phi|^2\rangle$ [$T_\mathrm{i}\rho_\mathrm{ti}/(eR_a)$]")
    ax.plot(kx,np.sqrt(ds_phikx.phi[it,:]),".-")
    ax.set_xlim(-4,4)
    ax.set_ylim(0,None)
    # ax.set_title(labels[it]+" t={:5.2f}".format(tt[it]))
    ax.set_title(labels[it])

fig.tight_layout()
plt.savefig("fig4_rotating.pdf",dpi=600,bbox_inches="tight")
plt.show()


# In[ ]:




