## Initialize
import Pkg
Pkg.activate("/global/homes/a/aramreye/.julia/environments/v1.5/")
using Plots, BenchmarkTools, Profile
gr()
Pkg.activate("/global/homes/a/aramreye/RamirezReyes_Yang_2020_SpontaneousCyclogenesis/Project.toml")
using AvailablePotentialEnergyFramework
using Statistics:mean
using NCDatasets
import ImageSegmentation: seeded_region_growing, labels_map, region_adjacency_graph

## Read test file
test_file_with_cyclones = "/global/cscratch1/sd/aramreye/sam3d/subsetsForLin/f5e-4_2km_1000km_control_2d.nc"
test_file_with_cyclones_3d = "/global/cscratch1/sd/aramreye/sam3d/subsetsForLin/f5e-4_2km_1000km_control_3d.nc"
ds = Dataset(test_file_with_cyclones);
ds_3d = Dataset(test_file_with_cyclones_3d);
pres = ds["PSFC"][:,:,end- 10:end] :: Array{Float32,3};
U = ds_3d["U"][:,:,:,end - 10:end] :: Array{Float32,4};
V = ds_3d["V"][:,:,:, end - 10: end] :: Array{Float32,4};
x = ds_3d["x"][:] :: Array{Float32,1};
z = ds_3d["z"][:] :: Array{Float32,1};
close(ds)
close(ds_3d) 
windspeed = hypot.(U,V);
pres_anomaly = pres[:,:,:] .- mean(pres[:,:,:],dims=(1,2));

function speedtest1()
  test_file_with_cyclones_3d = "/global/cscratch1/sd/aramreye/sam3d/subsetsForLin/f5e-4_2km_1000km_control_3d.nc"  
  ds_3d = Dataset(test_file_with_cyclones_3d);
  U = ds_3d["U"][:,:,:,end - 10:end] :: Array{Float32,4};
  close(ds_3d)   
end

function speedtest2()
  test_file_with_cyclones_3d = "/global/cscratch1/sd/aramreye/sam3d/subsetsForLin/f5e-4_2km_1000km_control_3d.nc"  
  ds_3d = Dataset(test_file_with_cyclones_3d);
  U = variable(ds_3d,"U")[:,:,:,end - 10:end] :: Array{Float32,4};
  close(ds_3d)   
end


?variable

@btime speedtest1()

@btime speedtest2()

## Test center detection
function test_centerdetection(pres_anomaly,windspeed)
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
end



test_centerdetection(pres_anomaly,windspeed)

## Test cyclone detection and segmentation
function test_cyclonesegmentation(pressure,pres_anomaly)
    anim = @animate for t in 1:size(windspeed,4)
        centers_and_labels,cyclones = AvailablePotentialEnergyFramework.detect_cyclones(pres[:,:,t],-5,2000)
        labels_map2 = copy(labels_map(cyclones));
        labels_map2[labels_map2 .== 1000] .= 0;
        p1 = @views heatmap(permutedims(pres_anomaly[:,:,t]),c=:viridis,title = "Pressure anomaly (hPa)")
        p2 = heatmap(permutedims(labels_map2),reuse=false,title="Segmented cyclones timestep = $t")
        plot(p1,p2,size=(800,380))
    end
    gif(anim, "testresults/cyclone_segmentation.gif", fps = 1)
end

test_cyclonesegmentation(pres,pres_anomaly)

#Test cyclone averaging
function test_cycloneaddition_withmask!(presaddition,windspeedaddition,buf2d,buf3d,pres,windspeed,pres_anomaly)
    totalcyclonecount = 0 
    anim = @animate for t in 1:size(windspeed,4)
        centers_and_labels,cyclones = AvailablePotentialEnergyFramework.detect_cyclones(pres[:,:,t],-5,2000)
        if !isnothing(centers_and_labels[1])
            count2d = AvailablePotentialEnergyFramework.add_allcyclones!(presaddition,buf2d,pres[:,:,t],cyclones,centers_and_labels)
            count3d = AvailablePotentialEnergyFramework.add_allcyclones!(windspeedaddition,buf3d,windspeed[:,:,:,t],cyclones,centers_and_labels)
            totalcyclonecount += count2d
            p1 = @views heatmap(1e-3x,1e-3x,permutedims(pres_anomaly[:,:,t]),c=:viridis,title = "Pressure anomaly (hPa)")
            p2 = @views heatmap(1e-3x,1e-3x,permutedims(presaddition[:,:])./totalcyclonecount,title="Averaged surface pressure (hPa)",clims=(950,1010))
            p3 = @views contourf(1e-3x,1e-3z,permutedims(windspeedaddition[:,256,:])./totalcyclonecount,title="Averaged wind speed (m/s)")
            plot(p1,p2,p3,size=(1200,390), layout = (1,3))
        end
    end
    gif(anim, "testresults/cycloneaveraging_withmask.gif", fps = 1)
end

size_pres = size(pres)[1:2]
size_windspeed = size(windspeed)[1:3]
presaddition = zeros(size_pres)
windspeedaddition = zeros(size_windspeed)
buf2d = similar(presaddition)
buf3d = similar(windspeedaddition)
test_cycloneaddition_withmask!(presaddition,windspeedaddition,buf2d,buf3d,pres,windspeed,pres_anomaly)

#Test cyclone averaging
function test_cycloneaddition_withoutmask!(presaddition,windspeedaddition,buf2d,buf3d,pres,windspeed,pres_anomaly)
    totalcyclonecount = 0 
    anim = @animate for t in 1:size(windspeed,4)
        centers_and_labels,cyclones = AvailablePotentialEnergyFramework.detect_cyclones(pres[:,:,t],-5,2000)
        if !isnothing(centers_and_labels[1])
            count2d = AvailablePotentialEnergyFramework.add_allcyclones!(presaddition,buf2d,pres[:,:,t],cyclones,centers_and_labels,false)
            count3d = AvailablePotentialEnergyFramework.add_allcyclones!(windspeedaddition,buf3d,windspeed[:,:,:,t],cyclones,centers_and_labels,false)
            totalcyclonecount += count2d
            p1 = @views heatmap(1e-3x,1e-3x,permutedims(pres_anomaly[:,:,t]),c=:viridis,title = "Pressure anomaly (hPa)")
            p2 = @views heatmap(1e-3x,1e-3x,permutedims(presaddition[:,:])./totalcyclonecount,title="Averaged surface pressure (hPa)",clims=(950,1010))
            p3 = @views contourf(1e-3x,1e-3z,permutedims(windspeedaddition[:,256,:])./totalcyclonecount,title="Averaged wind speed (m/s)")
            plot(p1,p2,p3,size=(1200,390), layout = (1,3))
        end
    end
    gif(anim, "testresults/cycloneaveraging_withoutmask.gif", fps = 1)
end

size_pres = size(pres)[1:2]
size_windspeed = size(windspeed)[1:3]
presaddition = zeros(size_pres)
windspeedaddition = zeros(size_windspeed)
buf2d = similar(presaddition)
buf3d = similar(windspeedaddition)
test_cycloneaddition_withoutmask!(presaddition,windspeedaddition,buf2d,buf3d,pres,windspeed,pres_anomaly)

#Test cyclone averaging
function test_radialaverage_withoutmask!(pressureradial,windspeedradial,presaddition,windspeedaddition,buf2d,buf3d,radiusbins,pres,windspeed,pres_anomaly)
    totalcyclonecount = 0 
    anim = @animate for t in 1:size(windspeed,4)
        centers_and_labels,cyclones = AvailablePotentialEnergyFramework.detect_cyclones(pres[:,:,t],-5,2000)
        if !isnothing(centers_and_labels[1])
            count2d = AvailablePotentialEnergyFramework.add_allcyclones!(presaddition,buf2d,pres[:,:,t],cyclones,centers_and_labels,false)
            count3d = AvailablePotentialEnergyFramework.add_allcyclones!(windspeedaddition,buf3d,windspeed[:,:,:,t],cyclones,centers_and_labels,false)
            totalcyclonecount += count2d
            for idx in 1:(length(radiusbins) - 1)
                windspeedradial[idx,:] .= AvailablePotentialEnergyFramework.averageallindistance( (radiusbins[idx],radiusbins[idx+1]),windspeedaddition./totalcyclonecount,(256,256), 2000 )
                pressureradial[idx] = AvailablePotentialEnergyFramework.averageallindistance( (radiusbins[idx],radiusbins[idx+1]),presaddition./totalcyclonecount,(256,256), 2000 )
            end
            p1 = @views heatmap(1e-3x,1e-3x,permutedims(pres_anomaly[:,:,t]),c=:viridis,title = "Pressure anomaly (hPa)")
            p2 = @views plot(1e-3(radiusbins[1:end-1]),pressureradial,title="Averaged surface pressure (hPa)",clims=(950,1010))
            p3 = @views contourf(1e-3(radiusbins[1:end-1]),1e-3z,permutedims(windspeedradial),title="Averaged wind speed (m/s)")
            plot(p1,p2,p3,size=(1200,390), layout = (1,3))
        end
    end
    gif(anim, "testresults/composite_radial_withoutmask.gif", fps = 1)
end

#Test azimuthal averaging
size_pres = size(pres)[1:2]
size_windspeed = size(windspeed)[1:3]
presaddition = zeros(size_pres)
windspeedaddition = zeros(size_windspeed)
buf2d = similar(presaddition)
buf3d = similar(windspeedaddition)
radiusbins = 1000:2000:512000;
presradial = zeros(length(radiusbins) - 1)
windspeedradial = zeros(length(radiusbins) - 1,size(windspeed,3))
test_radialaverage_withoutmask!(presradial,windspeedradial,presaddition,windspeedaddition,buf2d,buf3d,radiusbins,pres,windspeed,pres_anomaly)

#Test cyclone averaging
function test_radialaverage_withmask!(pressureradial,windspeedradial,presaddition,windspeedaddition,buf2d,buf3d,radiusbins,pres,windspeed,pres_anomaly)
    totalcyclonecount = 0 
    anim = @animate for t in 1:size(windspeed,4)
    #for t in 1:size(windspeed,4)
        centers_and_labels,cyclones = AvailablePotentialEnergyFramework.detect_cyclones(pres[:,:,t],-5,2000)
        if !isnothing(centers_and_labels[1])
            #@info centers_and_labels[1]
            count2d = AvailablePotentialEnergyFramework.add_allcyclones!(presaddition,buf2d,pres[:,:,t],cyclones,centers_and_labels,true)
            count3d = AvailablePotentialEnergyFramework.add_allcyclones!(windspeedaddition,buf3d,windspeed[:,:,:,t],cyclones,centers_and_labels,true)
            totalcyclonecount += count2d
            for idx in 1:(length(radiusbins) - 1)
                if !iszero(totalcyclonecount)
                    windspeedradial[idx,:] .= AvailablePotentialEnergyFramework.averageallindistance( (radiusbins[idx],radiusbins[idx+1]),windspeedaddition./totalcyclonecount,(256,256), 2000 )
                    pressureradial[idx] = AvailablePotentialEnergyFramework.averageallindistance( (radiusbins[idx],radiusbins[idx+1]),presaddition./totalcyclonecount,(256,256), 2000 )
                end
            end
            p1 = @views heatmap(1e-3x,1e-3x,permutedims(pres_anomaly[:,:,t]),c=:viridis,title = "Pressure anomaly (hPa)")
            p2 = @views plot(1e-3.*(radiusbins[1:end-1]),pressureradial,title="Averaged surface pressure (hPa)",clims=(950,1010))
            p3 = @views contourf(1e-3.*(radiusbins[1:end-1]),1e-3z,permutedims(windspeedradial),title="Averaged wind speed (m/s)")
            plot(p1,p2,p3,size=(1200,390), layout = (1,3))
        end
    end
    gif(anim, "testresults/composite_radial_withmask.gif", fps = 1)
end

#Test azimuthal averaging
size_pres = size(pres)[1:2]
size_windspeed = size(windspeed)[1:3]
presaddition = zeros(size_pres)
windspeedaddition = zeros(size_windspeed)
buf2d = similar(presaddition)
buf3d = similar(windspeedaddition)
radiusbins = 1000:2000:512000;
presradial = ones(length(radiusbins) - 1)
windspeedradial = ones(length(radiusbins) - 1,size(windspeed,3))
test_radialaverage_withmask!(presradial,windspeedradial,presaddition,windspeedaddition,buf2d,buf3d,radiusbins,pres,windspeed,pres_anomaly)


