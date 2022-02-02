using AvailablePotentialEnergyFramework, PyCall,PyPlot, Statistics, JLD; tkr = pyimport("matplotlib.ticker")# import FormatStrFormatter

#using , Statistics, DataStructures
#using JLD, SAMtools, ProgressMeter

#const ApeHelperFunctions = SAMtools

data_dir = "/global/cscratch1/sd/aramreye/for_postprocessing/ApeBudgetOutputs/"

# apes_1daysmooth_f5e-4_2km_1000km_homoAll_homoSfc_part_1.jld
function gatherapes(data_dir,exp_name)
            int_ape = Float32[]
            int_ape_ke = Float32[]
            int_ape_dia = Float32[]
            int_ape_rad = Float32[]
            int_ape_ub2 =Float32[]
            int_ape_vb2 = Float32[]
            int_ape_wn2 = Float32[]
            int_ape_rate = Float32[]
            int_ape_Fs = Float32[]
    for iter in 1:24
    file_name = string(exp_name,iter,".jld")
    file_path = joinpath(data_dir,file_name)

    
        file = jldopen(file_path, "r", mmaparrays=true)
            append!(int_ape,read(file, "int_APE"))
            append!(int_ape_ke,read(file, "int_KE"))
            append!(int_ape_dia,read(file, "int_DIA"))
            append!(int_ape_rad,read(file, "int_RAD"))
            append!(int_ape_ub2,read(file, "int_Ub2"))
            append!(int_ape_vb2,read(file, "int_Vb2"))
            append!(int_ape_wn2,read(file, "int_WN2"))
            append!(int_ape_rate,read(file, "int_APE_rate"))
            append!(int_ape_Fs,read(file, "APE_Fs"))
        close(file)
    end
    residual = int_ape_rate .-  int_ape_rad .-  int_ape_dia .+  int_ape_ub2 .+ int_ape_vb2 .+  int_ape_wn2 .- int_ape_Fs
    return int_ape, int_ape_ke, int_ape_rad, int_ape_dia, int_ape_ub2.+int_ape_vb2, int_ape_wn2, int_ape_Fs, residual, int_ape_rate
end


alldata = [gatherapes(data_dir,"f5e-4_2km_1000km_control")
 gatherapes(data_dir,"f5e-4_2km_1000km_homoRad")
gatherapes(data_dir,"f5e-4_2km_1000km_homoSfc")
gatherapes(data_dir,"f5e-4_2km_1000km_homoRad_homoSfc")];

int_APE_control = gatherapes(data_dir,"f5e-4_2km_1000km_control")
int_APE_homoRad = gatherapes(data_dir,"f5e-4_2km_1000km_homoRad")
int_APE_homoSfc = gatherapes(data_dir,"f5e-4_2km_1000km_homoSfc")
int_APE_homoRad_homoSfc = gatherapes(data_dir,"f5e-4_2km_1000km_homoRad_homoSfc");

size(int_APE_homoRad_homoSfc)

smooth_window = 12
smooth_window_2 = 10
int_APE_control_smooth = filter_array_time.(int_APE_control,smooth_window,1)
int_APE_homoRad_smooth = filter_array_time.(int_APE_homoRad,smooth_window,1)
int_APE_homoSfc_smooth = filter_array_time.(int_APE_homoSfc,smooth_window,1) 
int_APE_homoRad_homoSfc_smooth = filter_array_time.(int_APE_homoRad_homoSfc,smooth_window,1);



times = range(0,stop=100,length=1200)

experiment_labels = ["Control" "HomoSfc" "HomoRad" "HomoAll"]
budget_labels = ["APE" "KE" "Rad" "Conv" "Advec" "WN2" "Fs" "Residual" "APE change rate"]
plot_labels = ["a" "b" "c" "d" "e"]

PyPlot.matplotlib.rc("font", family="sans-serif",size=14)
legend_fontsize, plotlabel_fontsize = 12,16
#legend_fontsize, plotlabel_fontsize = 6,8

fig, ax = plt.subplots(5,1,sharex=true,figsize=(18,9))
ax[1].plot(times,1e-5int_APE_control_smooth[1],label = experiment_labels[1])
ax[1].plot(times,1e-5int_APE_homoSfc_smooth[1],label = experiment_labels[2])
ax[1].plot(times,1e-5int_APE_homoRad_smooth[1],label = experiment_labels[3])
ax[1].plot(times,1e-5int_APE_homoRad_homoSfc_smooth[1],label = experiment_labels[4])
ax[1].set_ylabel(string("APE ", L"(J/m^2)",L"x10^5"))
ax[1].yaxis.set_label_coords(-0.05,0.5)
ax[1].yaxis.set_major_formatter(tkr.FormatStrFormatter("%d"))

for variable in 3:7
   ax[2].plot(times,int_APE_control_smooth[variable],label=budget_labels[variable])
end    

for variable in 3:7
   ax[3].plot(times,int_APE_homoSfc_smooth[variable])
end 
for variable in 3:7
   ax[4].plot(times,int_APE_homoRad_smooth[variable])
end 
for variable in 3:7
   ax[5].plot(times,int_APE_homoRad_homoSfc_smooth[variable])
end 
for plotid in 1:5
    ax[plotid].annotate(plot_labels[plotid], xy=(0.02, 0.8), xycoords="axes fraction",backgroundcolor="white")
    ax[plotid].grid(b=true,which="both",color="xkcd:gray", linestyle="--",alpha=0.2)
    if plotid > 1
        ax[plotid].text(-0.0, 1.03, experiment_labels[plotid-1], transform=ax[plotid,1].transAxes, size=legend_fontsize)
    end
end

ax[3].set_ylabel(string("APE budget ", L"(W/m^2)"))
ax[3].yaxis.set_label_coords(-0.05,0.25)      
ax[5].set_xlabel("Time (days)")


ax[1].legend(loc = (0.1, 1), ncol=4 ,frameon=false,fontsize=legend_fontsize) 
ax[2].legend(loc = (0.1, 1), ncol=8 ,frameon=false,fontsize=legend_fontsize)
ax[1].xaxis.set_minor_locator(tkr.MultipleLocator(5))
ax[1].yaxis.set_minor_locator(plt.NullLocator())


