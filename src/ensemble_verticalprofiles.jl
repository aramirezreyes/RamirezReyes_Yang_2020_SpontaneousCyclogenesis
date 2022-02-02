using Markdown
using InteractiveUtils

using Pkg
Pkg.activate("/global/u2/a/aramreye/RamirezReyes_Yang_2020_SpontaneousCyclogenesis")
using PyPlot
using PyCall
using Base64
using Printf
using Statistics
using AvailablePotentialEnergyFramework
using DelimitedFiles

file_control = "/global/homes/a/aramreye/SAM_model/SAM6.10.6/exp/f5e-4_2km_1000km_control/init_snd/snd300"
file_ensemble1 = "/global/homes/a/aramreye/SAM_model/SAM6.10.6/exp/f5e-4_2km_1000km_control_ensemble1/init_snd/snd300"
file_ensemble2 = "/global/homes/a/aramreye/SAM_model/SAM6.10.6/exp/f5e-4_2km_1000km_control_ensemble2/init_snd/snd300"

profile_control = readdlm(file_control;header=true,skipstart=1)[1][1:56,2:4]
profile_ensemble1= readdlm(file_ensemble1;header=true,skipstart=1)[1][1:56,2:4]
profile_ensemble2 = readdlm(file_ensemble2;header=true,skipstart=1)[1][1:56,2:4]

profile_enseble1_modified = copy(profile_ensemble1)

PyPlot.matplotlib.rc("font", family="sans-serif",size=8)

fig, ax = plt.subplots(1,2,sharey=true,figsize=(3,4),dpi=400)

ax[1].plot(profile_control[:,2],profile_control[:,1],ls)
ax[1].plot(profile_ensemble1[:,2],profile_ensemble1[:,1])
ax[1].plot(profile_ensemble2[:,2],profile_ensemble2[:,1])
ax[1].invert_yaxis()
ax[1].set_ylabel(string("Pressure ",L"(hPa)"))
ax[1].set_xlabel(string("Potential Temperature ",L"(K)"))
ax[2].plot(profile_control[:,3],profile_control[:,1])
ax[2].plot(profile_ensemble1[:,3],profile_ensemble1[:,1])
ax[2].plot(profile_ensemble2[:,3],profile_ensemble2[:,1])
ax[2].set_ylabel(string("Specific humidiy ",L"(g/kg)"))
plt.tight_layout()


