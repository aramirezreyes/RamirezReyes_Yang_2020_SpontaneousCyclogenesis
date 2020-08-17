module RamirezReyes_Yang_SpontaneousCyclogenesis

using ImageSegmentation
using AvailablePotentialEnergyFramework
using JLD
using NCDatasets
using Statistics



# Write your package code here.

include("APE_calculation_from_sam_output.jl")
include("Composite_creation_from_sam_output_and_diagnostics.jl")
export getapeanalysis,get_composites

#import animate_tests.jl

end
