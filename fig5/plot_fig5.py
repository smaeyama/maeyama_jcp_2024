#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

plt.style.use('../nature_style.txt')

ds_phirz = xr.open_dataset("data/rotating_phi_rzt.nc")
tt = np.array(ds_phirz.t)
rr = np.array(ds_phirz.rr)
zz= np.array(ds_phirz.zz)
phi = np.array(ds_phirz.phi)
# labels = ["(a)","(b)","(c)","(d)","(e)","(f)"]
# labels = [r"(a) $\bar{t}=0$",r"(b) $\bar{t}=T_\mathrm{lap}/8$",r"(c) $\bar{t}=T_\mathrm{lap}/4$",r"(d) $\bar{t}=T_\mathrm{lap}/2$",r"(e) $\bar{t}=3T_\mathrm{lap}/4$",r"(d) $\bar{t}=T_\mathrm{lap}$"]
labels = [r"(a) $t=T_\mathrm{lap}$",r"(b) $t=9T_\mathrm{lap}/8$",r"(c) $t=5T_\mathrm{lap}/4$",r"(d) $t=3T_\mathrm{lap}/2$",r"(e) $t=7T_\mathrm{lap}/4$",r"(d) $t=2T_\mathrm{lap}$"] 

fig=plt.figure(figsize=(7,4.0),dpi=600) # figsize=(width,height(inch)),dpi(dots per inch)
it=0
ax=fig.add_subplot(2,3,it+1)
ax.set_ylabel(r"Height $Z/R_\mathrm{a}$")
val = 0.6/np.exp(0.131*tt[0])
vmax = val*np.exp(0.131*tt[it])
quad = ax.pcolormesh(rr[it,:,:],zz[it,:,:],phi[it,:,:],cmap="jet",vmax=vmax,vmin=-vmax, rasterized=True)
# ax.set_title(labels[it]+" t={:5.2f}".format(tt[it]))
ax.set_title(labels[it])
fig.colorbar(quad,shrink=0.89)
ax.set_aspect("equal")
for it in range(1,3):
    ax=fig.add_subplot(2,3,it+1)
    vmax = val*np.exp(0.131*tt[it])
    quad = ax.pcolormesh(rr[it,:,:],zz[it,:,:],phi[it,:,:],cmap="jet",vmax=vmax,vmin=-vmax, rasterized=True)
    # ax.set_title(labels[it]+" t={:5.2f}".format(tt[it]))
    ax.set_title(labels[it])
    fig.colorbar(quad,shrink=0.89)
    ax.set_aspect("equal")
it=3
ax=fig.add_subplot(2,3,it+1)
ax.set_xlabel(r"Major radius $R/R_\mathrm{a}$")
ax.set_ylabel(r"Height $Z/R_\mathrm{a}$")
val = 0.6/np.exp(0.131*tt[0])
vmax = val*np.exp(0.131*tt[it])
quad = ax.pcolormesh(rr[it,:,:],zz[it,:,:],phi[it,:,:],cmap="jet",vmax=vmax,vmin=-vmax, rasterized=True)
# ax.set_title(labels[it]+" t={:5.2f}".format(tt[it]))
ax.set_title(labels[it])
fig.colorbar(quad,shrink=0.89)
ax.set_aspect("equal")
for it in range(4,6):
    ax=fig.add_subplot(2,3,it+1)
    ax.set_xlabel(r"Major radius $R/R_\mathrm{a}$")
    vmax = val*np.exp(0.131*tt[it])
    quad = ax.pcolormesh(rr[it,:,:],zz[it,:,:],phi[it,:,:],cmap="jet",vmax=vmax,vmin=-vmax, rasterized=True)
    # ax.set_title(labels[it]+" t={:5.2f}".format(tt[it]))
    ax.set_title(labels[it])
    fig.colorbar(quad,shrink=0.89)
    ax.set_aspect("equal")


fig.tight_layout()
plt.savefig("fig5.pdf",dpi=600,bbox_inches="tight")
plt.show()


# In[ ]:




