#!/bin/bash -l
#SBATCH --qos=premium
#SBATCH --time=15:00:00
#SBATCH --nodes=3
#SBATCH --constraint=haswell
#SBATCH --mail-user=aramirezreyes@ucdavis.edu
#SBATCH --license=project,SCRATCH
#SBATCH --mail-type=begin
#SBATCH --mail-type=end



export TMPDIR=$SCRATCH

#/global/homes/a/aramreye/Software/julia-1.5.1/bin/julia --project=@. -e 'using RamirezReyes_Yang_SpontaneousCyclogenesis; using Distributed; using ClusterManagers; addprocs(SlurmManager(9),m="cyclic",exeflags="--project"); @everywhere using RamirezReyes_Yang_SpontaneousCyclogenesis; RamirezReyes_Yang_SpontaneousCyclogenesis.computebudgets("f5e-4_2km_1000km_control")' 


/global/homes/a/aramreye/Software/julia-1.5.1/bin/julia --project=@. -e 'using RamirezReyes_Yang_SpontaneousCyclogenesis; using Distributed; using ClusterManagers; addprocs(SlurmManager(9),m="cyclic",exeflags="--project=/global/u2/a/aramreye/RamirezReyes_Yang_2020_SpontaneousCyclogenesis"); @everywhere using RamirezReyes_Yang_SpontaneousCyclogenesis; RamirezReyes_Yang_SpontaneousCyclogenesis.computebudgets_last15days("f5e-4_2km_1000km_control")'

/global/homes/a/aramreye/Software/julia-1.5.1/bin/julia --project=@. -e 'using RamirezReyes_Yang_SpontaneousCyclogenesis; using Distributed; using ClusterManagers; addprocs(SlurmManager(9),m="cyclic",exeflags="--project=/global/u2/a/aramreye/RamirezReyes_Yang_2020_SpontaneousCyclogenesis"); @everywhere using RamirezReyes_Yang_SpontaneousCyclogenesis; RamirezReyes_Yang_SpontaneousCyclogenesis.computebudgets_last15days("f5e-4_2km_1000km_homoRad")'

/global/homes/a/aramreye/Software/julia-1.5.1/bin/julia --project=@. -e 'using RamirezReyes_Yang_SpontaneousCyclogenesis; using Distributed; using ClusterManagers; addprocs(SlurmManager(9),m="cyclic",exeflags="--project=/global/u2/a/aramreye/RamirezReyes_Yang_2020_SpontaneousCyclogenesis"); @everywhere using RamirezReyes_Yang_SpontaneousCyclogenesis; RamirezReyes_Yang_SpontaneousCyclogenesis.computebudgets_last15days("f5e-4_2km_1000km_homoSfc")'

/global/homes/a/aramreye/Software/julia-1.5.1/bin/julia --project=@. -e 'using RamirezReyes_Yang_SpontaneousCyclogenesis; using Distributed; using ClusterManagers; addprocs(SlurmManager(9),m="cyclic",exeflags="--project=/global/u2/a/aramreye/RamirezReyes_Yang_2020_SpontaneousCyclogenesis"); @everywhere using RamirezReyes_Yang_SpontaneousCyclogenesis; RamirezReyes_Yang_SpontaneousCyclogenesis.computebudgets_last15days("f5e-4_2km_1000km_homoRad_homoSfc")'
