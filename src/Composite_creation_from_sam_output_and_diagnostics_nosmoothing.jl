

"""
        getcomposites(file2d,file3d,diagnostics_dir,exp_,initial_timeindex,final_timeindex,partnumber,outputInterval)


"""
function get_composites_nosmoothing(file2d,file3d,exp_name)

    @info "Starting composite routine"
    diagnostics_data_dir = "/global/cscratch1/sd/aramreye/for_postprocessing/ApeBudgetOutputs_nosmoothing/"
    composites_data_dir = "/global/cscratch1/sd/aramreye/for_postprocessing/CompositeOutputs_nosmoothing/"
    azimuthal_averages_data_dir = "/global/cscratch1/sd/aramreye/for_postprocessing/AzimuthualAverageOutputs_nosmoothing/"
    
    variables_diag = ["convec_heating_anomaly","rad_heating_anomaly","buoyancy_anomaly"]
    variables_orig_3d = ["U","V","QV","TABS","QRAD","PP","W"]
    variables_orig_2d = ["USFC","VSFC","PW","Prec","LHF","SHF"]
    
    initial_timeindex = 601
    #    initial_timeindex = 1195
    final_timeindex = 1200
    
    iterator_time_2d    = initial_timeindex*2-1:2:final_timeindex*2
    iterator_time_3d    = initial_timeindex:1:final_timeindex
    
    radiusbins = 1000:2000:512000; 
    ds3d                = Dataset(file3d)
    ds2d                = Dataset(file2d)
    z                   = variable(ds3d,"z")[:]                       :: Array{Float32,1}
    t                   = variable(ds3d,"time")[iterator_time_3d]     :: Array{Float32,1}
    surface_pressure    = variable(ds2d,"PSFC")[:,:,iterator_time_2d] :: Array{Float32,3}
    pressure_anomaly    = surface_pressure .- mean(surface_pressure,dims=(1,2))

    ##Create buffers to reuse on allocating functions

    size_var_of_interest_2d = size(surface_pressure)[1:end-1]
    size_var_of_interest_3d = (size_var_of_interest_2d...,length(z))
    ## Create 3D buffers
    #var_of_interest_3d      = zeros(size_var_of_interest_3d)
    composite_addition_3d   = zeros(Float32,size_var_of_interest_3d)
    buf_3d_2                = similar(composite_addition_3d)
    buf_3d_1                = similar(composite_addition_3d)
    azimuthal_average_3d    = zeros(Float32,length(radiusbins) - 1,length(z))
    ## Create 2D buffers
    #var_of_interest_2d      = zeros(size_var_of_interest_2d)
    composite_addition_2d   = zeros(Float32,size_var_of_interest_2d)
    buf_2d_1                = similar(composite_addition_2d)
    buf_2d_2                = similar(composite_addition_2d)
    azimuthal_average_2d    = zeros(Float32,length(radiusbins) - 1)
    @info size_var_of_interest_2d size_var_of_interest_3d
    


    centers_labels_and_cyclones = @views [
        detect_cyclones!(buf_2d_1, buf_2d_2, pressure_anomaly[:,:,timeindex],-5,2000) for
        timeindex in 1:size(surface_pressure,3)]; pressure_anomaly = [];

    composite_filename_withmask = string(composites_data_dir,exp_name        ,"_withmask.jld")
    composite_filename_nomask   = string(composites_data_dir,exp_name        ,"_nomask.jld")
    azimuthal_filename_withmask = string(azimuthal_averages_data_dir,exp_name,"_withmask.jld")
    azimuthal_filename_nomask   = string(azimuthal_averages_data_dir,exp_name,"_nomask.jld")


    ### *** add allcyclones ***
    @info "Creating composite for surface pressure "
    flush(stdout)
    composite_addition_2d .= 0.0
    buf_2d_1              .= 0.0
    buf_2d_2              .= 0.0
    azimuthal_average_2d  .= 0.0
    composite_and_radialaverage!(azimuthal_average_2d,composite_addition_2d,buf_2d_1,
    buf_2d_2,centers_labels_and_cyclones,radiusbins,surface_pressure;maskcyclones = true)

    @info "Writing"
    flush(stdout)
    h5open(composite_filename_withmask, "w") do file
        write(file,"surface_pressure",composite_addition_2d)
    end
    

    h5open(azimuthal_filename_withmask, "w") do file
        write(file,"surface_pressure",azimuthal_average_2d)
    end
    GC.gc()
    GC.gc()
    ## Create products without mask
    @info "Creating composite without mask for: surface_pressure"
    flush(stdout)
    composite_addition_2d .= 0.0
    buf_2d_1              .= 0.0
    buf_2d_2              .= 0.0
    azimuthal_average_2d  .= 0.0
    composite_and_radialaverage!(azimuthal_average_2d,composite_addition_2d,buf_2d_1,buf_2d_2
    ,centers_labels_and_cyclones,radiusbins,surface_pressure;maskcyclones = false)
    @info "Writing"
    flush(stdout)
    h5open(composite_filename_nomask, "w") do file
        write(file,"surface_pressure",composite_addition_2d)
    end        
    h5open(azimuthal_filename_nomask, "w") do file
        write(file,"surface_pressure",azimuthal_average_2d)
    end
    GC.gc()
    GC.gc()
    

    

    for current_variable in variables_orig_3d
        ##NCDatasets.load!(ds["PSFC"].var,a,:,:,:)
        
        var_of_interest = variable(ds3d,current_variable)[:,:,:,iterator_time_3d]    :: Array{Float32,4}
        composite_addition_3d .= 0.0
        buf_3d_1              .= 0.0
        buf_3d_2              .= 0.0
        azimuthal_average_3d  .= 0.0
        

        ## Create products with mask
        @info "Creating composite for: " current_variable
        flush(stdout)
        composite_and_radialaverage!(azimuthal_average_3d,composite_addition_3d,buf_3d_1,buf_3d_2,centers_labels_and_cyclones,radiusbins,var_of_interest;maskcyclones = true)
        @info "Writing"
        flush(stdout)
        h5open(composite_filename_withmask, "r+",mmaparrays=true) do file
            write(file,current_variable,composite_addition_3d)
        end        
        h5open(azimuthal_filename_withmask, "r+",mmaparrays=true) do file
            write(file,current_variable,azimuthal_average_3d)
        end
        GC.gc()
        GC.gc()
        ## Create products without mask
        @info "Creating composite without mask for: " current_variable
        flush(stdout)
        composite_addition_3d .= 0.0
        buf_3d_1              .= 0.0
        buf_3d_2              .= 0.0
        azimuthal_average_3d  .= 0.0
        composite_and_radialaverage!(azimuthal_average_3d,composite_addition_3d,buf_3d_1,buf_3d_2,centers_labels_and_cyclones,radiusbins,var_of_interest;maskcyclones = false)
        @info "Writing"
        flush(stdout)
        h5open(composite_filename_nomask, "r+",mmaparrays=true) do file
            write(file,current_variable,composite_addition_3d)
        end        
        h5open(azimuthal_filename_nomask, "r+",mmaparrays=true) do file
            write(file,current_variable,azimuthal_average_3d)
        end
        GC.gc()
    end
    

    for current_variable in variables_orig_2d
        composite_addition_2d .= 0.0
        buf_2d_1              .= 0.0
        buf_2d_2              .= 0.0
        azimuthal_average_2d  .= 0.0

        var_of_interest = variable(ds2d,current_variable)[:,:,iterator_time_2d]    :: Array{Float32,3}    
        
        ## Create products with mask
        @info "Creating composite for: " current_variable
        flush(stdout)
        composite_and_radialaverage!(azimuthal_average_2d,composite_addition_2d,buf_2d_1,buf_2d_2,centers_labels_and_cyclones,radiusbins,var_of_interest;maskcyclones = true)
        @info "Writing"
        flush(stdout)
        h5open(composite_filename_withmask, "r+",mmaparrays=true) do file
            write(file,current_variable,composite_addition_2d)
        end        
        h5open(azimuthal_filename_withmask, "r+",mmaparrays=true) do file
            write(file,current_variable,azimuthal_average_2d)
        end
        GC.gc()
        GC.gc()
        ## Create products without mask
        @info "Creating composite without mask for: " current_variable
        flush(stdout)
        composite_addition_2d .= 0.0
        buf_2d_1              .= 0.0
        buf_2d_2              .= 0.0
        azimuthal_average_2d  .= 0.0
        composite_and_radialaverage!(azimuthal_average_2d,composite_addition_2d,buf_2d_1,buf_2d_2,centers_labels_and_cyclones,radiusbins,var_of_interest;maskcyclones = false)
        @info "Writing"
        flush(stdout)
        h5open(composite_filename_nomask, "r+",mmaparrays=true) do file
            write(file,current_variable,composite_addition_2d)
        end        
        h5open(azimuthal_filename_nomask, "r+",mmaparrays=true) do file
            write(file,current_variable,azimuthal_average_2d)
        end
        GC.gc()
    end


    for current_variable in variables_diag
        var_of_interest = gather_diagnostic_3d(diagnostics_data_dir,exp_name,current_variable)   :: Array{Float32,4}
        composite_addition_3d .= 0.0
        buf_3d_1              .= 0.0
        buf_3d_2              .= 0.0
        azimuthal_average_3d  .= 0.0
        ## Create products with mask
        @info "Creating composite for: " current_variable
        flush(stdout)
        composite_and_radialaverage!(azimuthal_average_3d,composite_addition_3d,buf_3d_1,buf_3d_2,centers_labels_and_cyclones,radiusbins,var_of_interest;maskcyclones = true)
        @info "Writing"
        flush(stdout)
        h5open(composite_filename_withmask, "r+",mmaparrays=true) do file
            write(file,current_variable,composite_addition_3d)
        end        
        h5open(azimuthal_filename_withmask, "r+",mmaparrays=true) do file
            write(file,current_variable,azimuthal_average_3d)
        end
        GC.gc()
        GC.gc()
        ## Create products without mask
        @info "Creating composite without mask for: " current_variable
        flush(stdout)
        composite_addition_3d .= 0.0
        buf_3d_1              .= 0.0
        buf_3d_2              .= 0.0
        azimuthal_average_3d  .= 0.0
        composite_and_radialaverage!(azimuthal_average_3d,composite_addition_3d,buf_3d_1,buf_3d_2,centers_labels_and_cyclones,radiusbins,var_of_interest;maskcyclones = false)
        @info "Writing"
        flush(stdout)
        h5open(composite_filename_nomask, "r+",mmaparrays=true) do file
            write(file,current_variable,composite_addition_3d)
        end        
        h5open(azimuthal_filename_nomask, "r+",mmaparrays=true) do file
            write(file,current_variable,azimuthal_average_3d)
        end
        GC.gc()
        GC.gc()
    end

    

    close(ds3d)   
    close(ds2d)   

    
    return nothing

end


"""
        getcomposites(file2d,file3d,diagnostics_dir,exp_,initial_timeindex,final_timeindex,partnumber,outputInterval)


"""
function cyclone_detect_and_save(file2d,file3d,exp_name)

    @info "Starting composite routine"
    diagnostics_data_dir = "/global/cscratch1/sd/aramreye/for_postprocessing/ApeBudgetOutputs_nosmoothing/"
    composites_data_dir = "/global/cscratch1/sd/aramreye/for_postprocessing/CompositeOutputs_nosmoothing/"
    azimuthal_averages_data_dir = "/global/cscratch1/sd/aramreye/for_postprocessing/AzimuthualAverageOutputs_nosmoothing/"
    
    variables_diag = ["convec_heating_anomaly","rad_heating_anomaly","buoyancy_anomaly"]
    variables_orig_3d = ["U","V","QV","TABS","QRAD","PP","W"]
    variables_orig_2d = ["USFC","VSFC","PW","Prec","LHF","SHF"]
    
    initial_timeindex = 601
    #    initial_timeindex = 1195
    final_timeindex = 620
#    final_timeindex = 1200
    
    iterator_time_2d    = initial_timeindex*2-1:2:final_timeindex*2
    iterator_time_3d    = initial_timeindex:1:final_timeindex
    
    radiusbins = 1000:2000:512000; 
    ds3d                = Dataset(file3d)
    ds2d                = Dataset(file2d)
    z                   = variable(ds3d,"z")[:]                       :: Array{Float32,1}
    t                   = variable(ds3d,"time")[iterator_time_3d]     :: Array{Float32,1}
    surface_pressure    = variable(ds2d,"PSFC")[:,:,iterator_time_2d] :: Array{Float32,3}
    pressure_anomaly    = surface_pressure .- mean(surface_pressure,dims=(1,2))

    ##Create buffers to reuse on allocating functions

    size_var_of_interest_2d = size(surface_pressure)[1:end-1]
    size_var_of_interest_3d = (size_var_of_interest_2d...,length(z))
    ## Create 3D buffers
    #var_of_interest_3d      = zeros(size_var_of_interest_3d)
    composite_addition_3d   = zeros(Float32,size_var_of_interest_3d)
    buf_3d_2                = similar(composite_addition_3d)
    buf_3d_1                = similar(composite_addition_3d)
    azimuthal_average_3d    = zeros(Float32,length(radiusbins) - 1,length(z))
    ## Create 2D buffers
    #var_of_interest_2d      = zeros(size_var_of_interest_2d)
    composite_addition_2d   = zeros(Float32,size_var_of_interest_2d)
    buf_2d_1                = similar(composite_addition_2d)
    buf_2d_2                = similar(composite_addition_2d)
    azimuthal_average_2d    = zeros(Float32,length(radiusbins) - 1)
    @info size_var_of_interest_2d size_var_of_interest_3d
    

    for timeindex in 1:size(surface_pressure,3)
        centers_labels_and_cyclones = @views detect_cyclones!(buf_2d_1, buf_2d_2, pressure_anomaly[:,:,timeindex],-5,2000)
        @info "Data types " typeof(centers_labels_and_cyclones[1]) centers_labels_and_cyclones[1]  typeof(centers_labels_and_cyclones[2])
    end
        pressure_anomaly = [];

#    jldopen(string(composites_data_dir,"centers_detection_",exp_name,".jld"), "w") do file
#        addrequire(file, AvailablePotentialEnergyFramework)
#        write(file, "centers_data", centers_labels_and_cyclones)  # alternatively, say "@write file A"
#    end

    

    composite_filename_withmask = string(composites_data_dir,exp_name        ,"_withmask.jld")
    composite_filename_nomask   = string(composites_data_dir,exp_name        ,"_nomask.jld")
    azimuthal_filename_withmask = string(azimuthal_averages_data_dir,exp_name,"_withmask.jld")
    azimuthal_filename_nomask   = string(azimuthal_averages_data_dir,exp_name,"_nomask.jld")

    
    

    close(ds3d)   
    close(ds2d)   

    
    return nothing

end
