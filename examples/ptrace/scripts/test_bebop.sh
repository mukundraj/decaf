#!/bin/bash

# executable
exe=./ptrace


# max_step
dtSim=$((5000*7200*10))

# stepsize
dtParticle=300000

# seed rate
seed_rate=100

# sample percentage
pred_percent=10

# skip rate / skip percentage
skip_rate=10

# without prediction
prediction=0


# num_procs=4
# particle_file="particles/particles4.nc"
# gp_file="graphinfo/graph.info.part.4"
# opts="--blocks 4 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle"
# mpiexec -n $num_procs $exe $args

# num_procs=8
# particle_file="particles/particles8.nc"
# gp_file="graphinfo/graph.info.part.8"
# opts="--blocks 8 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle"
# mpiexec -n $num_procs $exe $args

num_procs=16
particle_file="particles/particles16.nc"
gp_file="graphinfo/graph.info.part.16"
opts="--blocks 16 --check 1"
args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent $skip_rate"
mpiexec -n $num_procs $exe $args

# num_procs=128
# particle_file="particles/particle128.nc"
# gp_file="graphinfo/graph.info.part.128"
# opts="--blocks 128 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle"
# mpiexec -n $num_procs $exe $args



# with prediction
prediction=1


# num_procs=4
# particle_file="particles/particles4.nc"
# gp_file="graphinfo/graph.info.part.4"
# opts="--blocks 4 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle"
# mpiexec -n $num_procs $exe $args

# num_procs=8
# particle_file="particles/particles8.nc"
# gp_file="graphinfo/graph.info.part.8"
# opts="--blocks 8 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle"
# mpiexec -n $num_procs $exe $args


num_procs=16
particle_file="particles/particles16.nc"
gp_file="graphinfo/graph.info.part.16"
opts="--blocks 16 --check 1"
args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent $skip_rate"
mpiexec -n $num_procs $exe $args

