#!/bin/bash -l
#SBATCH --clusters=escori
#SBATCH --qos=bigmem
#SBATCH --nodes=1
#SBATCH --time=14:00:00
#SBATCH --mail-user=aramirezreyes@ucdavis.edu
#SBATCH --license=project,SCRATCH
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mem=300GB

export JULIA_NUM_THREADS=1
export TMPDIR=$SCRATCH

/global/homes/a/aramreye/Software/julia-1.5.0/bin/julia --project=@. -e 'using RamirezReyes_Yang_SpontaneousCyclogenesis; smooth_vars_and_write_to_netcdf!("/global/cscratch1/sd/aramreye/for_postprocessing/largencfiles/smoothed_variables/f5e-4_2km_1000km_homoRad_3d_smoothed.nc","/global/cscratch1/sd/aramreye/for_postprocessing/largencfiles/f5e-4_2km_1000km_homoRad_3d.nc",("U","V", "W", "QV", "TABS", "QRAD","PP"),11,60)' 

#/global/homes/a/aramreye/Software/julia-1.5.0/bin/julia --project=@. -e 'using RamirezReyes_Yang_SpontaneousCyclogenesis; smooth_vars_and_write_to_netcdf!("/global/cscratch1/sd/aramreye/sam3d/subsetsForLin/f5e-4_2km_1000km_control_3d_smoothed.nc","/global/cscratch1/sd/aramreye/sam3d/subsetsForLin/f5e-4_2km_1000km_control_3d.nc",("U","V", "W", "QV", "TABS", "QRAD"),11,60)' 



#("U","V", "W", "QV", "TABS", "QRAD")

