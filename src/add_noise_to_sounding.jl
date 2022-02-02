using SAMTools
using AvailablePotentialEnergyFramework
using DelimitedFiles

file_control = "/global/homes/a/aramreye/SAM_model/SAM6.10.6/exp/f5e-4_2km_1000km_control/init_snd/snd300"
file_ensemble1 = "/global/homes/a/aramreye/SAM_model/SAM6.10.6/exp/f5e-4_2km_1000km_control_ensemble1/init_snd/snd300"
file_ensemble2 = "/global/homes/a/aramreye/SAM_model/SAM6.10.6/exp/f5e-4_2km_1000km_control_ensemble2/init_snd/snd300"

profile_control = readdlm(file_control;header=true,skipstart=1)[1][1:56,2:4]
profile_ensemble1= readdlm(file_ensemble1;header=true,skipstart=1)[1][1:56,2:4]
profile_ensemble2 = readdlm(file_ensemble2;header=true,skipstart=1)[1][1:56,2:4]

profile_enseble1_modified = copy(profile_ensemble1)