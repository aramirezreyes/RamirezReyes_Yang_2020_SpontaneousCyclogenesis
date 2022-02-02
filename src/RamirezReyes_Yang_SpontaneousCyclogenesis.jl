module RamirezReyes_Yang_SpontaneousCyclogenesis

using Distributed: @everywhere, @spawnat, fetch
#using ClusterManagers
using ImageSegmentation
using AvailablePotentialEnergyFramework
using JLD
using HDF5
using NCDatasets
using Statistics

#addprocs(SlurmManager(9),m="cyclic")

# Write your package code here.

include("APE_calculation_from_sam_output.jl")
include("APE_calculation_from_sam_output_nosmoothing.jl")
include("APE_calculation_from_sam_output_addcondensateloading.jl")
include("Composite_creation_from_sam_output_and_diagnostics.jl")
include("Composite_creation_from_sam_output_and_diagnostics_nosmoothing.jl")
include("Composite_creation_from_sam_output_and_diagnostics_in_chunks_9hpa_15days.jl")

export getapeanalysis,getapeanalysis_nosmoothing,get_composites, getapeanalysis_last15days
export get_composites_nosmoothing
export smooth_vars_and_write_to_netcdf!, cyclone_detect_and_save, get_composites_in_chunks
#import animate_tests.jl

end
