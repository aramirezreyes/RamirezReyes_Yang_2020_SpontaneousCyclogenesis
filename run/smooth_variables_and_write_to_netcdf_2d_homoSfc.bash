#!/bin/bash -l
#SBATCH --qos=regular
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --constraint=knl
#SBATCH --mail-user=aramirezreyes@ucdavis.edu
#SBATCH --license=project,SCRATCH
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


export JULIA_NUM_THREADS=1
export TMPDIR=$SCRATCH

/global/homes/a/aramreye/Software/julia-1.5.0/bin/julia --project=@. -e 'using RamirezReyes_Yang_SpontaneousCyclogenesis; smooth_vars_and_write_to_netcdf!("/global/cscratch1/sd/aramreye/for_postprocessing/largencfiles/f5e-4_2km_1000km_homoSfc_2d_smoothed.nc","/global/cscratch1/sd/aramreye/for_postprocessing/largencfiles/f5e-4_2km_1000km_homoSfc_2d.nc",("SHF","LHF","PSFC"),11,120)' 



#("U","V", "W", "QV", "TABS", "QRAD")

