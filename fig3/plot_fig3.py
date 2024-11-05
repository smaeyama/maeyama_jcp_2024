#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt

s_hat=0.8
gamma_e=0.2
T_lap=2*np.pi*s_hat/gamma_e

plt.style.use('../nature_style.txt')

data_gamma0 = np.loadtxt("data/gamma0_gkvp.eng.001")
data_remap_mj07nx31 = np.loadtxt("data/remap_mj07nx31_gkvp.eng.001")
data_remap_mj01nx04 = np.loadtxt("data/remap_mj01nx04_gkvp.eng.001")
data_rotating_mj07nx31 = np.loadtxt("data/rotating_mj07nx31_gkvp.eng.001")
data_rotating_mj01nx04 = np.loadtxt("data/rotating_mj01nx04_gkvp.eng.001")

# plt.style.use('../nature_style.txt')
plt.rcParams["xtick.top"] = False

fig=plt.figure(figsize=(3.5,2.5),dpi=600) # figsize=(width,height(inch)),dpi(dots per inch)
ax=fig.add_subplot(111)
ax.set_xlabel(r"Time $t$ [$R_\mathrm{a}/v_\mathrm{ti}$]")
ax.set_ylabel(r"Potential $\langle|\tilde{\phi}|^2\rangle$ [$(T_\mathrm{i}\rho_\mathrm{ti}^2/(eR_a))^2$]")
ax.axvline(T_lap,c="gray",lw=0.3,linestyle="--")
ax.axvline(2*T_lap,c="gray",lw=0.3,linestyle="--")
ax.plot([3*T_lap,3*T_lap],[30,8e7],c="gray",lw=0.3,linestyle="--")# ax.axvline(3*T_lap,c="gray",lw=0.3,linestyle="--")
ax.plot([4*T_lap,4*T_lap],[30,8e7],c="gray",lw=0.3,linestyle="--")# ax.axvline(4*T_lap,c="gray",lw=0.3,linestyle="--")
ax.plot(data_gamma0[:,0], data_gamma0[:,1],linestyle=(0,(1,1)),label=r"$\gamma_E=0$")
ax.plot(data_remap_mj07nx31[:,0], data_remap_mj07nx31[:,1],"-",label=r"Remap ($k_{x,\mathrm{min}}=0.36$)")
ax.plot(data_remap_mj01nx04[:,0], data_remap_mj01nx04[:,1],linestyle=(0,(5,1,1,1)),label=r"Remap ($k_{x,\mathrm{min}}=2.51$)")
ax.plot(data_rotating_mj07nx31[:,0], data_rotating_mj07nx31[:,1],linestyle=(0,(5,1)),c="k",label=r"Rotating ($k_{x,\mathrm{min}}=0.36$)")
ax.plot(data_rotating_mj01nx04[:,0], data_rotating_mj01nx04[:,1],linestyle=(0,(5,1,1,1,1,1)),label=r"Rotating ($k_{x,\mathrm{min}}=2.51$)")
ax.set_yscale("log")
ax.set_xlim(0,110)
ax.set_ylim(1e-5,1e8)
# ax.legend(loc="upper left", bbox_to_anchor=(1.13,1.02),edgecolor="k")
ax.legend(loc="lower right")

def forward(time):
    return time
def inverse(time_normalized):
    return time_normalized
secax = ax.secondary_xaxis("top", functions=(forward, inverse))
# secax.set_xlabel(r"Shifted time $\bar{t}=t-T_\mathrm{lap}$")
secax.set_xticks(np.linspace(0,4*T_lap,5))
secax.set_xticklabels(["0",r"$T_\mathrm{lap}$",r"$2T_\mathrm{lap}$",r"$3T_\mathrm{lap}$",r"$4T_\mathrm{lap}$"])
secax.minorticks_off()


fig.tight_layout()
plt.savefig("fig3.pdf",dpi=600,bbox_inches="tight")
plt.show()


# In[ ]:





# In[ ]:




