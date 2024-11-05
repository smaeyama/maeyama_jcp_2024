#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt

### Read time-series data files
### Calculate time-average and its standard error by block bootstrap method
def block_bootstrap(data, block_size, n_iterations=1000, step_size=None):
    if step_size is None:
        step_size = block_size

    blocks = []
    for start in range(0, len(data) - block_size + 1, step_size):
        blocks.append(data[start:start + block_size])
    n_blocks = len(blocks)
    block_means = np.array([np.mean(block) for block in blocks])
    n_resample = min(n_blocks, int(n_blocks * step_size / block_size))
    bootstrap_means = np.empty(n_iterations)
    for i in range(n_iterations):
        sampled_block_means = np.random.choice(block_means, size=n_resample, replace=True)
        bootstrap_means[i] = np.mean(sampled_block_means)

    mean_value = np.mean(bootstrap_means)
    std_error = np.std(bootstrap_means)

    return mean_value, std_error

def calc_mean_and_std_error(data_list, block_size, step_size=None):
    mean_list = []
    std_error_list = []
    for data in data_list:
        mean, std_error = block_bootstrap(data, block_size, step_size=step_size)
        mean_list.append(mean)
        std_error_list.append(std_error)
    return np.array(mean_list), np.array(std_error_list)

def normalized_mean_and_std(ave,std):
    '''
    val = ave +- std
    normalized_val = val/val[0]
                   ~ ave/ave[0] + (std - std[0]*ave/ave[0])/ave[0] + O(std^2)
    '''
    normalized_mean = ave/ave[0]
    normalized_std = (std + std[0]*ave/ave[0])/ave[0]
    return normalized_mean, normalized_std


# In[ ]:


qes_gamma000 = np.loadtxt("data/qes_gamma0.00.dat")
qes_gamma005 = np.loadtxt("data/qes_gamma0.05.dat")
qes_gamma010 = np.loadtxt("data/qes_gamma0.10.dat")
qes_gamma020 = np.loadtxt("data/qes_gamma0.20.dat")
qes_gamma030 = np.loadtxt("data/qes_gamma0.30.dat")
qes_gamma035 = np.loadtxt("data/qes_gamma0.35.dat")
qes_gamma040 = np.loadtxt("data/qes_gamma0.40.dat")
qes_gamma050 = np.loadtxt("data/qes_gamma0.50.dat")
qes_gamma060 = np.loadtxt("data/qes_gamma0.60.dat")

growthrate = np.log(qes_gamma000[250,1]/qes_gamma000[240,1])/(2*(qes_gamma000[250,0]-qes_gamma000[240,0]))
print("Maximum linear growthrate:",growthrate)

# Time average
itsta=-5000
gamma = [0,0.05,0.1,0.2,0.3,0.35,0.4,0.5,0.6]
qes_list = [qes_gamma000[itsta:,1],
            qes_gamma005[itsta:,1],
            qes_gamma010[itsta:,1],
            qes_gamma020[itsta:,1],
            qes_gamma030[itsta:,1],
            qes_gamma035[itsta:,1],
            qes_gamma040[itsta:,1],
            qes_gamma050[itsta:,1],
            qes_gamma060[itsta:,1]]
block_size=300
qes_ave, qes_std = calc_mean_and_std_error(qes_list,block_size=block_size,step_size=block_size//2)
nqes_ave, nqes_std = normalized_mean_and_std(qes_ave, qes_std)


# In[ ]:


qes_remap_gamma000_nx055mj04 = np.loadtxt("data/qes_remap_gamma0.00_nx055mj04.dat")
qes_remap_gamma005_nx055mj04 = np.loadtxt("data/qes_remap_gamma0.05_nx055mj04.dat")
qes_remap_gamma010_nx055mj04 = np.loadtxt("data/qes_remap_gamma0.10_nx055mj04.dat")
qes_remap_gamma020_nx055mj04 = np.loadtxt("data/qes_remap_gamma0.20_nx055mj04.dat")
qes_remap_gamma030_nx055mj04 = np.loadtxt("data/qes_remap_gamma0.30_nx055mj04.dat")
qes_remap_gamma035_nx055mj04 = np.loadtxt("data/qes_remap_gamma0.35_nx055mj04.dat")
qes_remap_gamma040_nx055mj04 = np.loadtxt("data/qes_remap_gamma0.40_nx055mj04.dat")
qes_remap_gamma050_nx055mj04 = np.loadtxt("data/qes_remap_gamma0.50_nx055mj04.dat")
qes_remap_gamma060_nx055mj04 = np.loadtxt("data/qes_remap_gamma0.60_nx055mj04.dat")

# Time average
itsta=-5000
gamma_nx055mj04 = [0,0.05,0.1,0.2,0.3,0.35,0.4,0.5,0.6]
qes_remap_list_nx055mj04 = [qes_remap_gamma000_nx055mj04[itsta:,1],
                            qes_remap_gamma005_nx055mj04[itsta:,1],
                            qes_remap_gamma010_nx055mj04[itsta:,1],
                            qes_remap_gamma020_nx055mj04[itsta:,1],
                            qes_remap_gamma030_nx055mj04[itsta:,1],
                            qes_remap_gamma035_nx055mj04[itsta:,1],
                            qes_remap_gamma040_nx055mj04[itsta:,1],
                            qes_remap_gamma050_nx055mj04[itsta:,1],
                            qes_remap_gamma060_nx055mj04[itsta:,1]]
block_size=300
qes_remap_ave_nx055mj04, qes_remap_std_nx055mj04 = calc_mean_and_std_error(qes_remap_list_nx055mj04,block_size=block_size,step_size=block_size//2)
nqes_remap_ave_nx055mj04, nqes_remap_std_nx055mj04 = normalized_mean_and_std(qes_remap_ave_nx055mj04, qes_remap_std_nx055mj04)


# In[ ]:


qes_remap_gamma000_nx440mj32 = np.loadtxt("data/qes_remap_gamma0.00_nx440mj32.dat")
qes_remap_gamma005_nx440mj32 = np.loadtxt("data/qes_remap_gamma0.05_nx440mj32.dat")
qes_remap_gamma010_nx440mj32 = np.loadtxt("data/qes_remap_gamma0.10_nx440mj32.dat")
qes_remap_gamma020_nx440mj32 = np.loadtxt("data/qes_remap_gamma0.20_nx440mj32.dat")
qes_remap_gamma030_nx440mj32 = np.loadtxt("data/qes_remap_gamma0.30_nx440mj32.dat")
qes_remap_gamma040_nx440mj32 = np.loadtxt("data/qes_remap_gamma0.40_nx440mj32.dat")

# Time average
itsta=-2500
gamma_nx440mj32 = [0,0.05,0.1,0.2,0.3,0.4]
qes_remap_list_nx440mj32 = [qes_remap_gamma000_nx440mj32[itsta:,1],
                            qes_remap_gamma005_nx440mj32[itsta:,1],
                            qes_remap_gamma010_nx440mj32[itsta:,1],
                            qes_remap_gamma020_nx440mj32[itsta:,1],
                            qes_remap_gamma030_nx440mj32[itsta:,1],
                            qes_remap_gamma040_nx440mj32[itsta:,1]]
block_size=300
qes_remap_ave_nx440mj32, qes_remap_std_nx440mj32 = calc_mean_and_std_error(qes_remap_list_nx440mj32,block_size=block_size,step_size=block_size//2)
nqes_remap_ave_nx440mj32, nqes_remap_std_nx440mj32 = normalized_mean_and_std(qes_remap_ave_nx440mj32, qes_remap_std_nx440mj32)


# In[ ]:


qes_remap_gamma000_nx028mj02 = np.loadtxt("data/qes_remap_gamma0.00_nx028mj02.dat")
qes_remap_gamma000_nx042mj03 = np.loadtxt("data/qes_remap_gamma0.00_nx042mj03.dat")
qes_remap_gamma000_nx055mj04 = np.loadtxt("data/qes_remap_gamma0.00_nx055mj04.dat")
qes_remap_gamma000_nx110mj08 = np.loadtxt("data/qes_remap_gamma0.00_nx110mj08.dat")
qes_remap_gamma000_nx220mj16 = np.loadtxt("data/qes_remap_gamma0.00_nx220mj16.dat")
qes_remap_gamma000_nx440mj32 = np.loadtxt("data/qes_remap_gamma0.00_nx440mj32.dat")

# Time average
itsta=-5000
itsta2=-1500
Lx_list = [50,75,100,200,400,800]
qes_remap_list_gamma000 = [qes_remap_gamma000_nx028mj02[itsta:,1],
                           qes_remap_gamma000_nx042mj03[itsta:,1],
                           qes_remap_gamma000_nx055mj04[itsta:,1],
                           qes_remap_gamma000_nx110mj08[itsta:,1],
                           qes_remap_gamma000_nx220mj16[itsta2:,1],
                           qes_remap_gamma000_nx440mj32[itsta2:,1]]
block_size=300
qes_remap_ave_gamma000, qes_remap_std_gamma000 = calc_mean_and_std_error(qes_remap_list_gamma000,block_size=block_size,step_size=block_size//2)
nqes_remap_ave_gamma000, nqes_remap_std_gamma000 = normalized_mean_and_std(qes_remap_ave_gamma000, qes_remap_std_gamma000)


# In[ ]:


qes_remap_gamma005_nx028mj02 = np.loadtxt("data/qes_remap_gamma0.05_nx028mj02.dat")
qes_remap_gamma005_nx042mj03 = np.loadtxt("data/qes_remap_gamma0.05_nx042mj03.dat")
qes_remap_gamma005_nx055mj04 = np.loadtxt("data/qes_remap_gamma0.05_nx055mj04.dat")
qes_remap_gamma005_nx110mj08 = np.loadtxt("data/qes_remap_gamma0.05_nx110mj08.dat")
qes_remap_gamma005_nx220mj16 = np.loadtxt("data/qes_remap_gamma0.05_nx220mj16.dat")
qes_remap_gamma005_nx440mj32 = np.loadtxt("data/qes_remap_gamma0.05_nx440mj32.dat")

# Time average
itsta=-5000
itsta2=-1500
Lx_list = [50,75,100,200,400,800]
qes_remap_list_gamma005 = [qes_remap_gamma005_nx028mj02[itsta:,1],
                           qes_remap_gamma005_nx042mj03[itsta:,1],
                           qes_remap_gamma005_nx055mj04[itsta:,1],
                           qes_remap_gamma005_nx110mj08[itsta:,1],
                           qes_remap_gamma005_nx220mj16[itsta2:,1],
                           qes_remap_gamma005_nx440mj32[itsta2:,1]]
block_size=300
qes_remap_ave_gamma005, qes_remap_std_gamma005 = calc_mean_and_std_error(qes_remap_list_gamma005,block_size=block_size,step_size=block_size//2)
nqes_remap_ave_gamma005, nqes_remap_std_gamma005 = normalized_mean_and_std(qes_remap_ave_gamma005, qes_remap_std_gamma005)


# In[ ]:


qes_gamma005_nx028mj02 = np.loadtxt("data/qes_gamma0.05_nx028mj02.dat")
qes_gamma005_nx042mj03 = np.loadtxt("data/qes_gamma0.05_nx042mj03.dat")
qes_gamma005_nx055mj04 = np.loadtxt("data/qes_gamma0.05_nx055mj04.dat")
qes_gamma005_nx110mj08 = np.loadtxt("data/qes_gamma0.05_nx110mj08.dat")
qes_gamma005_nx220mj16 = np.loadtxt("data/qes_gamma0.05_nx220mj16.dat")
qes_gamma005_nx440mj32 = np.loadtxt("data/qes_gamma0.05_nx440mj32.dat")

# Time average
itsta=-5000
Lx_list = [50,75,100,200,400,800]
qes_list_gamma005 = [qes_gamma005_nx028mj02[itsta:,1],
                     qes_gamma005_nx042mj03[itsta:,1],
                     qes_gamma005_nx055mj04[itsta:,1],
                     qes_gamma005_nx110mj08[itsta:,1],
                     qes_gamma005_nx220mj16[itsta:,1],
                     qes_gamma005_nx440mj32[itsta:,1]]
block_size=300
qes_ave_gamma005, qes_std_gamma005 = calc_mean_and_std_error(qes_list_gamma005,block_size=block_size,step_size=block_size//2)
nqes_ave_gamma005, nqes_std_gamma005 = normalized_mean_and_std(qes_ave_gamma005, qes_std_gamma005)


# In[ ]:


def calc_spectrum(qes):
    ista=-6000
    t=qes[ista:,0]
    f_raw=qes[ista:,1]
    f=(f_raw-np.average(f_raw))*np.hanning(len(f_raw))
    # plt.plot(f_raw-np.average(f_raw));plt.plot(f);plt.show()
    dt=np.average(np.diff(t))
    nt=len(t)
    omega=2.0*np.pi*np.fft.rfftfreq(len(t),d=dt)
    fw=np.fft.rfft(f)
    return omega, np.abs(fw)

omega_gamma000, spectrum_gamma000 = calc_spectrum(qes_gamma000)
omega_gamma005, spectrum_gamma005 = calc_spectrum(qes_gamma005)
omega_gamma010, spectrum_gamma010 = calc_spectrum(qes_gamma010)
omega_gamma020, spectrum_gamma020 = calc_spectrum(qes_gamma020)
omega_gamma030, spectrum_gamma030 = calc_spectrum(qes_gamma030)
omega_gamma035, spectrum_gamma035 = calc_spectrum(qes_gamma035)
omega_gamma040, spectrum_gamma040 = calc_spectrum(qes_gamma040)
omega_gamma050, spectrum_gamma050 = calc_spectrum(qes_gamma050)
omega_remap_gamma000, spectrum_remap_gamma000 = calc_spectrum(qes_remap_gamma000_nx055mj04)
omega_remap_gamma005, spectrum_remap_gamma005 = calc_spectrum(qes_remap_gamma005_nx055mj04)
omega_remap_gamma010, spectrum_remap_gamma010 = calc_spectrum(qes_remap_gamma010_nx055mj04)
omega_remap_gamma020, spectrum_remap_gamma020 = calc_spectrum(qes_remap_gamma020_nx055mj04)
omega_remap_gamma030, spectrum_remap_gamma030 = calc_spectrum(qes_remap_gamma030_nx055mj04)
omega_remap_gamma035, spectrum_remap_gamma035 = calc_spectrum(qes_remap_gamma035_nx055mj04)
omega_remap_gamma040, spectrum_remap_gamma040 = calc_spectrum(qes_remap_gamma040_nx055mj04)
omega_remap_gamma050, spectrum_remap_gamma050 = calc_spectrum(qes_remap_gamma050_nx055mj04)


# In[ ]:


plt.style.use('../nature_style.txt')

fig=plt.figure(figsize=(3.35,3.9),dpi=600) # figsize=(width,height(inch)),dpi(dots per inch)
gs=fig.add_gridspec(nrows=2,ncols=1,height_ratios=(1,1))
ax=fig.add_subplot(gs[0,0])
ax.set_xlabel(r"Time $t$ [$R_\mathrm{a}/v_\mathrm{ti}$]")
ax.set_ylabel(r"Heat flux $Q_\mathrm{i}$ [$Q_\mathrm{gB}$]")
ax.plot(qes_gamma000[:,0], qes_gamma000[:,1],"-",label=r"$0$")
# ax.plot(qes_gamma005[:,0], qes_gamma005[:,1],"-",label=r"$0.05$")
ax.plot(qes_gamma010[:,0], qes_gamma010[:,1],linestyle=(0,(5,1)),label=r"$0.1$")#v_\mathrm{ti}/R_\mathrm{a}$")
# ax.plot(qes_gamma020[:,0], qes_gamma020[:,1],"-",label=r"$0.2$")
ax.plot(qes_gamma030[:,0], qes_gamma030[:,1],linestyle=(0,(1,1)),label=r"$0.3$")#v_\mathrm{ti}/R_\mathrm{a}$")
ax.plot(qes_gamma035[:,0], qes_gamma035[:,1],linestyle=(0,(5,1,1,1)),label=r"$0.35$")#v_\mathrm{ti}/R_\mathrm{a}$")
ax.plot(qes_gamma040[:,0], qes_gamma040[:,1],linestyle=(0,(5,1,1,1,1,1)),label=r"$0.4$")#v_\mathrm{ti}/R_\mathrm{a}$")
# ax.plot(qes_gamma050[:,0], qes_gamma050[:,1],"-",label=r"$0.5$")
# ax.plot(qes_gamma060[:,0], qes_gamma060[:,1],"-",label=r"$0.6$")
wx = np.linspace(28,38,200)
wy = np.exp(2*growthrate*wx)*0.05/np.exp(2*growthrate*25)
ax.plot(wx,wy,c="k",lw=0.5)
ax.set_xlim(0,1050)
ax.set_ylim(1e-1,1.2e2)
ax.set_yscale("log")
ax.axvline(qes_gamma005[0,0],c="k",lw=0.5,linestyle="--")
ax.legend(loc="lower right",bbox_to_anchor=(1,0.08))


ax=fig.add_subplot(gs[1,0])
ax.plot(omega_gamma010,spectrum_gamma010,"-",lw=0.8,label=r"Rotating ($\gamma_E=0.1v_\mathrm{ti}/R_a)$")
ax.plot(omega_gamma040,spectrum_gamma040,"--",lw=0.8,label=r"Rotating ($\gamma_E=0.4v_\mathrm{ti}/R_a$)")
ax.plot(omega_remap_gamma010,spectrum_remap_gamma010,":",lw=0.8,label=r"Remap ($\gamma_E=0.1v_\mathrm{ti}/R_a$)")
ax.plot(omega_remap_gamma040,spectrum_remap_gamma040,"-.",lw=0.8,label=r"Remap ($\gamma_E=0.4v_\mathrm{ti}/R_a$)")
s_hat=0.8; gamma_E=0.1; T_lap = 2*np.pi*s_hat/gamma_E
for i in np.arange(1,4):
    freq = i*2*np.pi/T_lap
    ax.axvline(freq,linestyle="-",ymin=0.92-np.log(i)*0.05,ymax=0.99-np.log(i)*0.05,lw=1.2,c="k")
s_hat=0.8; gamma_E=0.4; T_lap = 2*np.pi*s_hat/gamma_E
for i in np.arange(1,4):
    freq = i*2*np.pi/T_lap
    ax.axvline(freq,linestyle="-",ymin=0.5-np.log(i)*0.08,ymax=0.6-np.log(i)*0.08,lw=1.2,c="k")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"Frequency $\omega$ [$v_\mathrm{ti}/R_\mathrm{a}$]")
ax.set_ylabel(r"Spectrum of heat flux $|Q_{\mathrm{i},\omega}|$ [a.u]")
ax.set_xlim(omega_gamma040[2],omega_gamma040[-1])
ax.set_ylim(1e-6,None)
ax.legend()

fig.text(0.8,0.838,r"$\gamma_E R_a /v_\mathrm{ti}=$")
fig.text(0.89,0.93,"(a)",color="k",fontfamily="sans-serif",fontweight="bold",fontsize=8)
fig.text(0.89,0.445,"(b)",color="k",fontfamily="sans-serif",fontweight="bold",fontsize=8)

fig.tight_layout()
plt.savefig("fig8.pdf",dpi=600,bbox_inches="tight")
plt.show()


# In[ ]:


print(np.mean((qes_ave[0],qes_remap_ave_nx055mj04[0],qes_remap_ave_nx440mj32[0])))
plt.style.use('../nature_style.txt')

fig=plt.figure(figsize=(3.35,4.0),dpi=600) # figsize=(width,height(inch)),dpi(dots per inch)
gs=fig.add_gridspec(nrows=2,ncols=1,height_ratios=(1,1))
plt.rcParams["xtick.top"] = False

ax=fig.add_subplot(gs[0,0])
ax.set_xlabel(r"Normalized $E\times B$ shearing rate $\gamma_E/\gamma_\mathrm{max}$")
ax.set_ylabel(r"Heat flux $Q_\mathrm{i}(\gamma_E)/Q_\mathrm{i}(\gamma_E=0)$")
ax.errorbar(gamma/growthrate,nqes_ave,yerr=nqes_std,
            fmt="o",capsize=2,label=r"Rotating ($L_x=100\rho_\mathrm{ti}$)")
ax.errorbar(gamma_nx055mj04[1:]/growthrate,nqes_remap_ave_nx055mj04[1:],yerr=nqes_remap_std_nx055mj04[1:],
            fmt="^",capsize=2,label=r"Remap ($L_x=100\rho_\mathrm{ti}$)")
ax.errorbar(gamma_nx440mj32[1:]/growthrate,nqes_remap_ave_nx440mj32[1:],yerr=nqes_remap_std_nx440mj32[1:],
            fmt="x",capsize=2,label=r"Remap ($L_x=800\rho_\mathrm{ti}$)")
ax.set_xlim(0,4)
ax.set_ylim(0,1.2)
ax.set_yticks(0.2*np.arange(7))
ax.legend(loc="center right")
def forward(gamma_normalized):
    return gamma_normalized * growthrate
def inverse(gamma):
    return gamma / growthrate
secax = ax.secondary_xaxis("top", functions=(forward, inverse))
secax.set_xlabel(r"$E\times B$ shearing rate $\gamma_E$ [$v_\mathrm{ti}/R_a]$")
secax.set_xticks(np.linspace(0,0.8,9))
secax.set_xticklabels(["0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8"])
secax.minorticks_off()

plt.rcParams["xtick.top"] = False

ax=fig.add_subplot(gs[1,0])
ax.set_xlabel(r"Radial box size $L_x/\rho_\mathrm{ti}$")
ax.set_ylabel(r"Heat flux $Q_\mathrm{i}$ [$Q_\mathrm{gB}$]")
ax.errorbar(Lx_list,qes_remap_ave_gamma000,yerr=qes_remap_std_gamma000,fmt=":+",capsize=2,label=r"$\gamma_E=0$")
ax.errorbar(Lx_list,qes_ave_gamma005,yerr=qes_std_gamma005,fmt=":^",capsize=2,label=r"Rotating ($\gamma_E/\gamma_\mathrm{max}=0.24$)")
ax.errorbar(Lx_list,qes_remap_ave_gamma005,yerr=qes_remap_std_gamma005,fmt=":x",capsize=2,label=r"Remap ($\gamma_E/\gamma_\mathrm{max}=0.24$)")
ax.set_xlim(0.00001,850)
ax.set_ylim(0,37)
ax.legend()
# def forward(Lx):
#     ky_dominant = 0.2
#     s_hat = 0.8
#     gamma_e = 0.05
#     kxmin = 2*np.pi/Lx
#     T_remap = kxmin/(ky_dominant*gamma_e)
#     T_lap = 2*np.pi*s_hat/gamma_e
#     return T_remap/T_lap
# def inverse(normalized_T_remap):
#     ky_dominant = 0.2
#     s_hat = 0.8
#     gamma_e = 0.05
#     T_lap = 2*np.pi*s_hat/gamma_e
#     T_remap = normalized_T_remap * T_lap
#     kxmin = T_remap*(ky_dominant*gamma_e)
#     Lx = 2*np.pi/kxmin
#     # print(kxmin,T_remap,normalized_T_remap,T_lap,Lx)
#     return Lx
def forward(Lx):
    ky_dominant = 0.2
    s_hat = 0.8
    gamma_e = 0.05
    kxmin = 2*np.pi/Lx
    T_remap = kxmin/(ky_dominant*gamma_e)
    T_lap = 2*np.pi*s_hat/gamma_e
    return T_lap/T_remap
def inverse(normalized_T_lap):
    ky_dominant = 0.2
    s_hat = 0.8
    gamma_e = 0.05
    T_lap = 2*np.pi*s_hat/gamma_e
    Lx = 2*np.pi*normalized_T_lap/(T_lap*(ky_dominant*gamma_e))
    return Lx
secax = ax.secondary_xaxis("top", functions=(forward, inverse))
secax.set_xlabel(r"Remapping frequency $T_\mathrm{lap}/T_\mathrm{remap}(k_y\rho_\mathrm{ti}=0.2)$")
# secax.set_xticks(np.linspace(0,0.8,9))
# secax.set_xticklabels(["0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8"])
secax.minorticks_off()

fig.text(0.87,0.865,"(a)",color="k", fontfamily="sans-serif", fontweight="bold", fontsize=8)
fig.text(0.87,0.38,"(b)",color="k", fontfamily="sans-serif", fontweight="bold", fontsize=8)

fig.tight_layout()
plt.savefig("fig9.pdf",dpi=600,bbox_inches="tight")
plt.show()


# In[ ]:





# In[ ]:




