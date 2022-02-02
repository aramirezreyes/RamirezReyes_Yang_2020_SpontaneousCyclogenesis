"""
    getapeanalysis_nosmoothing(file2d,file3d,outputfile,initial_timeindex,final_timeindex,partnumber,outputInterval,position,FloatType::Type=Float64)
"""
function getapeanalysis_nosmoothing(file2d,file3d,outputfile,initial_timeindex,final_timeindex,partnumber,outputInterval,FloatType::Type=Float64)

    @info "Starting APE analysis routine"
    day           = 86400
    sst           = 300
    dt            = outputInterval
    
    dayLength     = 60*60*24÷outputInterval; #How many data points make one day
    #final_time    = floor(Int,final_time_days*dayLength) 
    #initial_time  = floor(Int,initial_time_days*dayLength)+1
    
    Pref          = 1000*1e2                 #Pa
    Ts            = sst                      #Sea surface temperature
    qs            = 25.7*1e-3 
    Tvs           = Ts*(1+epsilon*qs)
    c1            = (R/Dryair.cp)
    

    iterator_time_2d    = initial_timeindex*2-1:2:final_timeindex*2
    iterator_time_3d    = initial_timeindex:1:final_timeindex

    smooth_x      = smooth_y = 15 #it was 11
    smooth_time   = floor(Int,dayLength*5)+1 

    @info " Reading files: " file2d file3d
    
    ds3d                = Dataset(file3d)
    ds2d                = Dataset(file2d)
    x                   = variable(ds3d,"x")[:]    :: Array{Float32,1}
    y                   = variable(ds3d,"y")[:]    :: Array{Float32,1}
    z                   = variable(ds3d,"z")[:]    :: Array{Float32,1}
    t                   = variable(ds3d,"time")[iterator_time_3d] :: Array{Float32,1}
    P0                  = variable(ds3d,"p")[:]    :: Array{Float32,1}
    U                   = variable(ds3d,"U")[:,:,:,iterator_time_3d] :: Array{Float32,4}
    @info "Read 10%"    
    V                   = variable(ds3d,"V")[:,:,:,iterator_time_3d] :: Array{Float32,4}
    @info "Read 20%"    
    W                   = variable(ds3d,"W")[:,:,:,iterator_time_3d] :: Array{Float32,4}
    @info "Read 30%"    
    RAD                 = variable(ds3d,"QRAD")[:,:,:,iterator_time_3d] :: Array{Float32,4}
    @info "Read 40%"    
    T                   = variable(ds3d,"TABS")[:,:,:,iterator_time_3d] :: Array{Float32,4}
    @info "Read 50%"    
    Tv                  = variable(ds3d,"QV")[:,:,:,iterator_time_3d] :: Array{Float32,4}
    @info "Read 60%"    
    PP                  = variable(ds3d,"PP")[:,:,:,iterator_time_3d] :: Array{Float32,4}
    @info "Read 70%"                  
    SHF                 = variable(ds2d,"SHF")[:,:,iterator_time_2d] :: Array{Float32,3}
    @info "Read 80%"                 
    LHF                 = variable(ds2d,"LHF")[:,:,iterator_time_2d] :: Array{Float32,3}
    @info "Read 90%"

    close(ds3d)   
    close(ds2d)   
    


    #if 3 nodes:
    # node 1: workers 1, 2, 5,8
    # node 2: workers 3, 6, 9
    # node 3: workers 4, 7, 10
    

    ###### Processing
    
    @. SHF    = g/(1*Dryair.cp*Ts)*(SHF)                            
    @. LHF    = g/(1*Dryair.cp*Ts)*(epsilon*Dryair.cp*Ts/Liquidwater.Lv*LHF) 
    @. SHF    +=  LHF     # Now it is transformed

    
    ThetaV       = similar(T)
    xBar_Pt      = Array{eltype(T),4}(undef,1,1,size(PP,3),size(PP,4))
    xBar_Tv      = Array{eltype(T),4}(undef,1,1,size(Tv,3),size(Tv,4))
    xBar_ThetaV  = Array{eltype(T),4}(undef,1,1,size(ThetaV,3),size(ThetaV,4))
    #******
    dx           = x[2]-x[1]
    dy           = y[2]-y[1]                            # Grid size
    kz           = length(z)                            # vertical levels
    kx           = length(x)                            # # of horizonal grid points
    ky           = length(y)                            # # of horizonal grid points
    @.  RAD      = RAD/day                              # K/s #Heating rate per second
    @.  Tv       = (1 + 1e-3*epsilon*Tv)*T              # Virtual temperature
    @.  P0       = P0*1e2
    PP          .= PP .+ reshape(P0,(1,1,kz,1))
    ThetaV      .= Tv.*(P0[1]./PP).^c1 # Virtual potential temp
    mean!(xBar_Pt,PP)                                     
    mean!(xBar_Tv,Tv)              
    mean!(xBar_ThetaV,ThetaV)   
    var_Tv       =  Tv     .- xBar_Tv
    var_ThetaV   =  ThetaV .- xBar_ThetaV
    rho0         = dropdims(xBar_Pt./R./xBar_Tv,dims=(1,2))
    B            = g .* var_ThetaV./xBar_ThetaV
    @. RAD       = RAD*(g/xBar_Tv)                    # convert unit to buoyancy
    PP           = []
    
    N2           = compute_N2(xBar_Tv,z)

    LHF        = []
    T          = []
    ThetaV     = []
    Tv         = []
    var_ThetaV = []
    var_Tv     = []
    @info "Computing buoyancy budget"
    ### NOTE that SHF is now the sum, saving memory
    # Buoyancy budget
    dz          = 50
    @info size(B), size(RAD), size(SHF), size(U),size(V) ,size(W), size(N2), size(dx),size(dy), size(dz), size(dt), size(x),size(y), size(z), size(t)
    Diabatic_other = get_diabatic_as_residual_buoyancy(B, RAD, SHF, U,V ,W, N2, dx,dy, dz, dt)
    GC.gc()
    @info "Computing ape budget for the whole troposphere"

    # APE budget
    z_up        = 15000
    z_BL        = 2000

     (int_mass,
     int_KE,
     int_APE,
     int_APE_rate,
     int_APE_Ub2,
     int_APE_Vb2,
     int_APE_WN2,
     int_APE_RAD,
     int_APE_DIA,
     xBar_APE_Fs,
     residual) =  getapebudget(B, U,V, W, N2, RAD, SHF, Diabatic_other, rho0, x,y, z, t, dx,dy, dz, dt, z_up)


    Diabatic_other .= Diabatic_other .- mean(Diabatic_other,dims=(1,2)) #They are now perturbations
    RAD            .= RAD .- mean(RAD,dims=(1,2))
    B              .= B .- mean(B,dims=(1,2))



    dia_ape = Diabatic_other.*B
    rad_ape = RAD.*B   

    @info "WritingData"


    jldopen(string(outputfile,partnumber,".jld"), "w") do file
        write(file,"int_APE",int_APE)
        write(file,"int_KE",int_KE)
        write(file,"int_RAD",int_APE_RAD)
        write(file,"int_DIA",int_APE_DIA)
        write(file,"int_WN2",int_APE_WN2)
        write(file,"int_Ub2",int_APE_Ub2)
        write(file,"int_Vb2",int_APE_Vb2)
        write(file,"int_APE_rate",int_APE_rate)
        write(file,"APE_Fs",xBar_APE_Fs)
        write(file,"convec_heating_anomaly",Diabatic_other)
        write(file,"rad_heating_anomaly",RAD)
        write(file,"buoyancy_anomaly",B)
        write(file,"radiative_ape_production",dia_ape)
        write(file,"convective_ape_production",rad_ape)
    end
    GC.gc()
    return nothing

end

function computebudgets_nosmoothing(exp_name)

    #file2d = "/global/cscratch1/sd/aramreye/sam3d/subsetsForLin/f5e-4_2km_1000km_control_2d.nc"
    #file3d = "/global/cscratch1/sd/aramreye/sam3d/subsetsForLin/f5e-4_2km_1000km_control_3d.nc"
    outputfile = string("/global/cscratch1/sd/aramreye/for_postprocessing/ApeBudgetOutputs_nosmoothing/",exp_name,"_nosmoothing")
    input_folder = "/global/cscratch1/sd/aramreye/for_postprocessing/largencfiles/"
    file2d = string(input_folder,exp_name,"_2d.nc")
    file3d = string(input_folder,exp_name,"_3d.nc")
    total_days = 100
    output_timestep = 7200 #seconds
    indices_in_part = 50
    indices_in_day = 86400 ÷ output_timestep
    number_of_parts = total_days*indices_in_day÷indices_in_part

    for part_number in 1:number_of_parts
        initial_timestep = indices_in_part*(part_number-1) + 1
        final_timestep = indices_in_part*(part_number)
        @info "Starting part $part_number with indices from" initial_timestep final_timestep
        flush(stdout)
        flush(stderr)
        getapeanalysis_nosmoothing(file2d,file3d,outputfile,initial_timestep,final_timestep,part_number,output_timestep,Float32)

    end

end
