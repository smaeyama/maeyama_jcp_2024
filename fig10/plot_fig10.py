#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

plt.style.use('../nature_style.txt')

ds_phirz_gamma000 = xr.open_dataset("data/gamma0_phi_rzt.nc")
ds_phirz_gamma020 = xr.open_dataset("data/gamma0.2_phi_rzt.nc")
ds_phirz_gamma030 = xr.open_dataset("data/gamma0.3_phi_rzt.nc")
ds_phirz_gamma035 = xr.open_dataset("data/gamma0.35_phi_rzt.nc")
print(ds_phirz_gamma000)
print(ds_phirz_gamma020)
print(ds_phirz_gamma030)
print(ds_phirz_gamma035)


# In[ ]:


fig=plt.figure(figsize=(3.53,3.7),dpi=600) # figsize=(width,height(inch)),dpi(dots per inch)
ax=fig.add_subplot(2,2,1)
rr = np.array(ds_phirz_gamma000.rr)[0,:,:]
zz = np.array(ds_phirz_gamma000.zz)[0,:,:]
phi = np.array(ds_phirz_gamma000.phi)[0,:,:]
# vmax = max([np.max(phi),-np.min(phi)]) * 0.8
vmax=35
quad = ax.pcolormesh(rr,zz,phi,cmap="jet",vmax=vmax,vmin=-vmax, rasterized=True)
# ax.set_xlabel(r"Major radius $R/R_\mathrm{a}$")
ax.set_ylabel(r"Height $Z/R_\mathrm{a}$")
ax.set_title(r"(a) $\gamma_E=0$")
ax.set_aspect("equal")
ax.set_xticklabels([])
# fig.colorbar(quad,shrink=1)

ax=fig.add_subplot(2,2,2)
rr = np.array(ds_phirz_gamma020.rr)[0,:,:]
zz = np.array(ds_phirz_gamma020.zz)[0,:,:]
phi = np.array(ds_phirz_gamma020.phi)[0,:,:]
# vmax = max([np.max(phi),-np.min(phi)]) * 0.8
vmax=35
quad = ax.pcolormesh(rr,zz,phi,cmap="jet",vmax=vmax,vmin=-vmax, rasterized=True)
# ax.set_xlabel(r"Major radius $R/R_\mathrm{a}$")
# ax.set_ylabel(r"Height $Z/R_\mathrm{a}$")
ax.set_title(r"(b) $\gamma_E=0.2v_\mathrm{ti}/R_a$")
ax.set_aspect("equal")
ax.set_xticklabels([])
ax.set_yticklabels([])
# fig.colorbar(quad,shrink=1)

ax=fig.add_subplot(2,2,3)
rr = np.array(ds_phirz_gamma030.rr)[0,:,:]
zz = np.array(ds_phirz_gamma030.zz)[0,:,:]
phi = np.array(ds_phirz_gamma030.phi)[0,:,:]
# vmax = max([np.max(phi),-np.min(phi)]) * 0.8
vmax=35
quad = ax.pcolormesh(rr,zz,phi,cmap="jet",vmax=vmax,vmin=-vmax, rasterized=True)
ax.set_xlabel(r"Major radius $R/R_\mathrm{a}$")
ax.set_ylabel(r"Height $Z/R_\mathrm{a}$")
ax.set_title(r"(c) $\gamma_E=0.3v_\mathrm{ti}/R_a$")
ax.set_aspect("equal")
# fig.colorbar(quad,shrink=1)

ax=fig.add_subplot(2,2,4)
rr = np.array(ds_phirz_gamma035.rr)[0,:,:]
zz = np.array(ds_phirz_gamma035.zz)[0,:,:]
phi = np.array(ds_phirz_gamma035.phi)[0,:,:]
# vmax = max([np.max(phi),-np.min(phi)]) * 0.8
vmax=35
quad = ax.pcolormesh(rr,zz,phi,cmap="jet",vmax=vmax,vmin=-vmax, rasterized=True)
ax.set_xlabel(r"Major radius $R/R_\mathrm{a}$")
# ax.set_ylabel(r"Height $Z/R_\mathrm{a}$")
ax.set_title(r"(d) $\gamma_E=0.35v_\mathrm{ti}/R_a$")
ax.set_aspect("equal")
ax.set_yticklabels([])
# fig.colorbar(quad,shrink=1)
cax=fig.add_axes([0.265,-0.01,0.49,0.02])
cbar=fig.colorbar(quad,orientation="horizontal",cax=cax)
cbar.set_ticks(np.linspace(-30,30,7))

# fig.tight_layout()
plt.savefig("fig10.pdf",dpi=600,bbox_inches="tight")
plt.show()


# In[ ]:




