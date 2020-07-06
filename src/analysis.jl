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
pres = ds_3d["PP"][:,:,1,end] :: Array{Float32,2};
U = ds_3d["U"][:,:,:,end] :: Array{Float32,3};
V = ds_3d["V"][:,:,:,end] :: Array{Float32,3};
close(ds)
close(ds_3d) 
##Algorithm has to be: Detect cyclone, recenter, detect again and recenter. It has to be done twice to avoid centering some areas that are on the edge of a cyclone and are still detected as a local maximum

windspeed = sqrt.(U.*U .+ V.*V);
pres_anomaly = pres .- mean(pres,dims=(1,2));

heatmap(permutedims(windspeed[:,:,1]),title="Wind speed field (m/s)")
p1 = heatmap(permutedims(pres_anomaly))
centers = AvailablePotentialEnergyFramework.findcyclonecenters_aspressureminima(pres_anomaly,-5,2000);
centers_x= [centers[i][1] for i in 1:length(centers)];
centers_y = [centers[i][2] for i in 1:length(centers)];
scatter!(p1,centers_x,centers_y,markersize=10,markershape=:cross,c=:red)


centers_and_labels,cyclones = AvailablePotentialEnergyFramework.detect_cyclones(pres,-5,2000)
labels_map2 = copy(labels_map(cyclones));
labels_map2[labels_map2 .== 1000] .= 0;
p2 = heatmap(permutedims(labels_map2),reuse=false,title="Identified and segmented cyclones")

G, vert_map = region_adjacency_graph(cyclones, (i,j)->true);

radiusbins = 1000:2000:512000;

presaverages = []
for i in 1:(length(radiusbins) - 1 )
push!(presaverages,
AvailablePotentialEnergyFramework.averageallindistance(radiusbins[i:i+1],pres,labels_map2.==1,centers_and_labels[1][1],2000))
end

plot(radiusbins[1:end-1]./1e3,presaverages./100,title="Azimuthally averaged surface pressure anomaly (hPa)")

windspeedaverages = AvailablePotentialEnergyFramework.averageallindistance(radiusbins[1:2],windspeed,labels_map2.==2,centers_and_labels[2][1],2000)
for i in 2:(length(radiusbins) - 1 )
    global windspeedaverages
    windspeedaverages = hcat(windspeedaverages,
        AvailablePotentialEnergyFramework.averageallindistance(radiusbins[i:i+1],windspeed,labels_map2.==2,centers_and_labels[2][1],2000))
end
heatmap(windspeedaverages,title = "Azimuthally averaged wind speed for a random cyclone (m/s)")
size(windspeedaverages)


presaverage = AvailablePotentialEnergyFramework.azimuthalaverage_allcyclones(radiusbins,pres,cyclones,centers_and_labels,2000)
windspeedaverage = AvailablePotentialEnergyFramework.azimuthalaverage_allcyclones(radiusbins,windspeed,cyclones,centers_and_labels,2000)
plot(radiusbins[1:end-1]./1e3,1:80,presaverage./100,title="Azimuthally averaged surface pressure anomaly (hPa) 3 cyclones 1 timestep")
heatmap(radiusbins[1:end-1]./1e3,1:80,windspeedaverage,title="Azimuthally averaged surface windspeed (m/s) 3 cyclones 1 timestep")