"""

    composite_and_radialaverage_withmask(azimuthal_average,composite_addition,buf,radiusbins,surface_pressure,variable_of_interest)
"""


function composite_and_radialaverage!(azimuthal_average,composite_addition,buf,radiusbins,surface_pressure,variable_of_interest;withmask = true)
    totalcyclonecount = 0
    nd = ndims(variable_of_interest)
    for timeindex in 1:size(variable_of_interest,nd)
        centers_and_labels,cyclones = detect_cyclones(surface_pressure[:,:,timeindex],-5,2000)
        if !isnothing(centers_and_labels[1])
            count = add_allcyclones!(composite_addition,buf,selectdim(variable_of_interest,nd,timeindex),cyclones,centers_and_labels,withmask)
            totalcyclonecount += count
        end
    end
    if !iszero(totalcyclonecount)
        composite_addition ./= totalcyclonecount
        if nd == 4
            for idx in 1:(length(radiusbins) - 1)            
                azimuthal_average[idx,:] .= averageallindistance((radiusbins[idx],radiusbins[idx+1]),composite_addition,(256,256), 2000 )
            end
        elseif nd == 3
            for idx in 1:(length(radiusbins) - 1)            
                azimuthal_average[idx] = averageallindistance((radiusbins[idx],radiusbins[idx+1]),composite_addition,(256,256), 2000 )
            end
        end
        
    end
end



function gather_diagnostic_3d(data_dir,exp_name,desired_var)
    diagnostic = zeros(Float32,512,512,80,360)
    for iter in 1:8
        file_name = string(exp_name,iter+16,".jld")
        file_path = joinpath(data_dir,file_name)
        file = jldopen(file_path, "r", mmaparrays=true)
        index_initial = (iter-1)*(50)+1-40
        index_end = iter*50-40
        var_from_file = read(file, desired_var) 
         if iter == 1
             diagnostic[:,:,:,1:10] .=  var_from_file[:,:,:,end-9:end]
         else
             diagnostic[:,:,:,index_initial:index_end] .=  var_from_file[:,:,:,:]
         end
        close(file)
    end
    return diagnostic
end



"""
    getcomposites(file2d,file3d,diagnostics_dir,exp_,initial_timeindex,final_timeindex,partnumber,outputInterval)


"""
function get_composites(file2d,file3d,exp_name)

    @info "Starting composite routine"
    diagnostics_data_dir = "/global/cscratch1/sd/aramreye/for_postprocessing/ApeBudgetOutputs/"
    composites_data_dir = "/global/cscratch1/sd/aramreye/for_postprocessing/CompositeOutputs/"
    azimuthal_averages_data_dir = "/global/cscratch1/sd/aramreye/for_postprocessing/AzimuthualAverageOutputs/"
    
    variables_diag = ["convec_heating_anomaly","rad_heating_anomaly","buoyancy_anomaly"]
    variables_orig_3d = ["U","V","QV","TABS","QRAD","PP","W"]
    variables_orig_2d = ["USFC","VSFC","PW","Precip","LHF","SHF"]
    
    initial_timeindex = 841
#    initial_timeindex = 1195
    final_timeindex = 1200
    
    iterator_time_2d    = initial_timeindex*2-1:2:final_timeindex*2
    iterator_time_3d    = initial_timeindex:1:final_timeindex
    
    radiusbins = 1000:2000:512000; 
    ds3d                = Dataset(file3d)
    ds2d                = Dataset(file2d)
    z                   = variable(ds3d,"z")[:]    :: Array{Float32,1}
    t                   = variable(ds3d,"time")[iterator_time_3d] :: Array{Float32,1}
    surface_pressure           = variable(ds2d,"PSFC")[:,:,iterator_time_2d]    :: Array{Float32,3}
    size_var_of_interest_2d = size(surface_pressure)[1:end-1]
    size_var_of_interest_3d = (size_var_of_interest_2d...,length(z))
    @info size_var_of_interest_2d size_var_of_interest_3d
    
    composite_addition_3d = zeros(size_var_of_interest_3d)
    buf_3d = similar(composite_addition_3d)
    azimuthal_average_3d = zeros(length(radiusbins) - 1,length(z))

    composite_addition_2d = zeros(size_var_of_interest_2d)
    buf_2d = similar(composite_addition_2d)
    azimuthal_average_2d = zeros(length(radiusbins) - 1)

    
    composite_filename_withmask = string(composites_data_dir,exp_name,"_withmask.jld")
    composite_filename_nomask   = string(composites_data_dir,exp_name,"_nomask.jld")
    azimuthal_filename_withmask = string(azimuthal_averages_data_dir,exp_name,"_withmask.jld")
    azimuthal_filename_nomask   = string(azimuthal_averages_data_dir,exp_name,"_nomask.jld")


        @info "Creating composite for surface pressure " 
        composite_and_radialaverage!(azimuthal_average_2d,composite_addition_2d,buf_2d,radiusbins,surface_pressure,surface_pressure;withmask = true)
        @info "Writing"
        jldopen(composite_filename_withmask, "w") do file
            write(file,"surface_pressure",composite_addition_2d)
        end        
        jldopen(azimuthal_filename_withmask, "w") do file
            write(file,"surface_pressure",azimuthal_average_2d)
        end
        
        ## Create products without mask
        @info "Creating composite without mask for: surface_pressure" 
        composite_and_radialaverage!(azimuthal_average_2d,composite_addition_2d,buf_2d,radiusbins,surface_pressure,surface_pressure;withmask = false)
        @info "Writing"
        jldopen(composite_filename_nomask, "w") do file
            write(file,"surface_pressure",composite_addition_2d)
        end        
        jldopen(azimuthal_filename_nomask, "w") do file
            write(file,"surface_pressure",azimuthal_average_2d)
        end
    

    

    for current_variable in variables_orig_3d[1:1]
        var_of_interest = variable(ds3d,current_variable)[:,:,:,iterator_time_3d]    :: Array{Float32,4}
        composite_addition_3d .= 0.0
        buf_3d .= 0.0
        azimuthal_average_3d .= 0.0
        

        ## Create products with mask
        @info "Creating composite for: " current_variable
        composite_and_radialaverage!(azimuthal_average_3d,composite_addition_3d,buf_3d,radiusbins,surface_pressure,var_of_interest;withmask = true)
        @info "Writing"
        jldopen(composite_filename_withmask, "r+") do file
            write(file,current_variable,composite_addition_3d)
        end        
        jldopen(azimuthal_filename_withmask, "r+") do file
            write(file,current_variable,composite_addition_3d)
        end
        
        ## Create products without mask
        @info "Creating composite without mask for: " current_variable
        composite_and_radialaverage!(azimuthal_average_3d,composite_addition_3d,buf_3d,radiusbins,surface_pressure,var_of_interest;withmask = false)
        @info "Writing"
        jldopen(composite_filename_nomask, "r+") do file
            write(file,current_variable,composite_addition_3d)
        end        
        jldopen(azimuthal_filename_nomask, "r+") do file
            write(file,current_variable,azimuthal_average_3d)
        end
    end
    

    for current_variable in variables_orig_2d[1:1]
        composite_addition_2d .= 0.0
        buf_2d .= 0.0
        azimuthal_average_2d .= 0.0

        var_of_interest = variable(ds2d,current_variable)[:,:,iterator_time_2d]    :: Array{Float32,3}    
       
        ## Create products with mask
        @info "Creating composite for: " current_variable
        composite_and_radialaverage!(azimuthal_average_2d,composite_addition_2d,buf_2d,radiusbins,surface_pressure,var_of_interest;withmask = true)
        @info "Writing"
        jldopen(composite_filename_withmask, "r+") do file
            write(file,current_variable,composite_addition_2d)
        end        
        jldopen(azimuthal_filename_withmask, "r+") do file
            write(file,current_variable,azimuthal_average_2d)
        end
        
        ## Create products without mask
        @info "Creating composite without mask for: " current_variable
        composite_and_radialaverage!(azimuthal_average_2d,composite_addition_2d,buf_2d,radiusbins,surface_pressure,var_of_interest;withmask = false)
        @info "Writing"
        jldopen(composite_filename_nomask, "r+") do file
            write(file,current_variable,composite_addition_2d)
        end        
        jldopen(azimuthal_filename_nomask, "r+") do file
            write(file,current_variable,azimuthal_average_2d)
        end
    end


    for current_variable in variables_diag[1:1]
        var_of_interest = gather_diagnostic_3d(diagnostics_data_dir,exp_name,current_variable)   :: Array{Float32,4}
        composite_addition_3d .= 0.0
        buf_3d .= 0.0
        azimuthal_average_3d .= 0.0
        

        ## Create products with mask
        @info "Creating composite for: " current_variable
        composite_and_radialaverage!(azimuthal_average_3d,composite_addition_3d,buf_3d,radiusbins,surface_pressure,var_of_interest;withmask = true)
        @info "Writing"
        jldopen(composite_filename_withmask, "r+") do file
            write(file,current_variable,composite_addition_3d)
        end        
        jldopen(azimuthal_filename_withmask, "r+") do file
            write(file,current_variable,azimuthal_average_3d)
        end
        
        ## Create products without mask
        @info "Creating composite without mask for: " current_variable
        composite_and_radialaverage!(azimuthal_average_3d,composite_addition_3d,buf_3d,radiusbins,surface_pressure,var_of_interest;withmask = false)
        @info "Writing"
        jldopen(composite_filename_nomask, "r+") do file
            write(file,current_variable,composite_addition_3d)
        end        
        jldopen(azimuthal_filename_nomask, "r+") do file
            write(file,current_variable,azimuthal_average_3d)
        end
    end

  

    close(ds3d)   
    close(ds2d)   

    
    return nothing

end
