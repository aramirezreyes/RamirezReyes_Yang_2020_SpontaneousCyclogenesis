## Initialize
import Pkg
Pkg.activate("/Users/arreyes/.julia/environments/v1.5/")
using Plots, BenchmarkTools
gr()
Pkg.activate("/Users/arreyes/.julia/dev/RamirezReyes_Yang_SpontaneousCyclogenesis/")
using AvailablePotentialEnergyFramework
using Statistics:mean
using NCDatasets
import ImageSegmentation: seeded_region_growing, labels_map, region_adjacency_graph
## Read test file
test_file_with_cyclones = "/Users/arreyes/Documents/Research/Developement/testNetCDF/f5e-4_2km_1000km_homoRad_2d_timesteps.nc"
test_file_with_cyclones_3d = "/Users/arreyes/Documents/Research/Developement/testNetCDF/f5e-4_2km_1000km_homoRad_3d_timesteps.nc"
ds = Dataset(test_file_with_cyclones);
ds_3d = Dataset(test_file_with_cyclones_3d);
pres = ds_3d["PP"][:,:,:,:] :: Array{Float32,4};
U = ds_3d["U"][:,:,:,:] :: Array{Float32,4};
V = ds_3d["V"][:,:,:,:] :: Array{Float32,4};
close(ds)
close(ds_3d) 
windspeed = sqrt.(U.*U .+ V.*V);
pres_anomaly = pres[:,:,1,:] .- mean(pres[:,:,1,:],dims=(1,2));

## Test center detection
anim = @animate for t in 1:size(windspeed,4)
    centers = @views AvailablePotentialEnergyFramework.findcyclonecenters_aspressureminima(pres_anomaly[:,:,t],-5,2000);
    centers_x= [centers[i][1] for i in 1:length(centers)];
    centers_y = [centers[i][2] for i in 1:length(centers)];
    p1 = @views heatmap(permutedims(windspeed[:,:,1,t]),title="Wind speed (m/s)",c=:viridis)
    p2 = @views heatmap(permutedims(pres_anomaly[:,:,t]),c=:viridis,title = "Pressure anomaly (hPa) and centers , timestep = $t")
    scatter!(p2,centers_x,centers_y,markersize=10,markershape=:cross,c=:red)
    plot(p1,p2,size=(800,380))
end
gif(anim, "testresults/center_detection.gif", fps = 1)

## Test cyclone detection and segmentation

anim = @animate for t in 1:size(windspeed,4)
    centers_and_labels,cyclones = AvailablePotentialEnergyFramework.detect_cyclones(pres[:,:,1,t],-5,2000)
    labels_map2 = copy(labels_map(cyclones));
    labels_map2[labels_map2 .== 1000] .= 0;
    p1 = @views heatmap(permutedims(pres_anomaly[:,:,t]),c=:viridis,title = "Pressure anomaly (hPa)")
    p2 = heatmap(permutedims(labels_map2),reuse=false,title="Segmented cyclones timestep = $t")
    plot(p1,p2,size=(800,380))
end
gif(anim, "testresults/cyclone_segmentation.gif", fps = 1)



#Test azimuthal averaging

radiusbins = 1000:2000:512000;
anim = @animate for t in 1:size(windspeed,4)
    centers_and_labels,cyclones = AvailablePotentialEnergyFramework.detect_cyclones(pres[:,:,1,t],-5,2000)
    if t == 1
        global presaverage = AvailablePotentialEnergyFramework.azimuthalaverage_allcyclones(radiusbins,pres[:,:,1,t],cyclones,centers_and_labels,2000)
        global windspeedaverage = AvailablePotentialEnergyFramework.azimuthalaverage_allcyclones(radiusbins,windspeed[:,:,:,t],cyclones,centers_and_labels,2000)
    else
        presaverage = 0.5(presaverage .+ AvailablePotentialEnergyFramework.azimuthalaverage_allcyclones(radiusbins,pres[:,:,1,t],cyclones,centers_and_labels,2000))
        windspeedaverage = 0.5(windspeedaverage .+ AvailablePotentialEnergyFramework.azimuthalaverage_allcyclones(radiusbins,windspeed[:,:,:,t],cyclones,centers_and_labels,2000))
    end
    p1 = @views heatmap(permutedims(pres_anomaly[:,:,t]),c=:viridis,title = "Pressure anomaly (hPa)")
    p2 = plot(radiusbins[1:end-1]./1e3,presaverage./100,title="Azimuthally averaged surface pressure anomaly (hPa) ")
    p3 = heatmap(radiusbins[1:end-1]./1e3,1:80,windspeedaverage,title="Azimuthally averaged surface windspeed (m/s) t = $t")
    plot(p1,p2,p3,size=(1200,390))
end
gif(anim, "testresults/azimuthal_averages_alltimesteps.gif", fps = 1)