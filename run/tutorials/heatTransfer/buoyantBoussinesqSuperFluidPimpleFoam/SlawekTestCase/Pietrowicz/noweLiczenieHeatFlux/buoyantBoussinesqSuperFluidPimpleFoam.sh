#!/bin/bash

# nazwa zadania do 12 znakow
#PBS -N UnmUs
#PBS -q main
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=8:mpiprocs=8:mem=10000MB
#PBS -l software=openfoam-v1612
#PBS -m be

# wejscie do katalogu z ktorego zostalo uruchomione zadanie
cd $PBS_O_WORKDIR

module load openfoam/v1612+-gcc5.2.0
source /usr/local/OpenFOAM/OpenFOAM-v1612+/etc/bashrc

# compile turbulentHeatFluxTemperature
~/OpenFOAM/pblasiak-v1612+/src/TurbulenceModels/incompressible/turbulentTransportModel/turbulentHeatFluxTemperature/Allwmake > log.AllwmakeTurbulentHeatFluxTemperature

# compile Helium
~/OpenFOAM/pblasiak-v1612+/src/transportModels/incompressible/Allwmake > log.AllwmakeHelium

# compile buoyantBoussinesqSuperFluidPimpleFoam
~/OpenFOAM/pblasiak-v1612+/applications/solvers/heatTransfer/buoyantBoussinesqSuperFluidPimpleFoam/Allwmake > log.Allwmake

# run the simulation
./Allrun
