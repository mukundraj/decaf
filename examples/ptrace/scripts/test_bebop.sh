#!/bin/bash

# executable
exe=./ptrace


# max_step
dtSim=$((5000*7200))

# stepsize
dtParticle=$((300000/6))

# stepsize
dtParticle=300000

# seed rate
seed_rate=1


# with prediction
prediction=1

num_procs=512

pred_percent=5

skip_rate=1
particle_file="particles/particles16.nc"
gp_file="graphinfo/graph.info.part.16"
opts="--blocks 16 --check 0"
args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent $skip_rate"
mpiexec -n $num_procs $exe $args

