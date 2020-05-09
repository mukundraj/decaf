#!/bin/bash

#SBATCH --job-name=p1
#SBATCH --account=pedal
#SBATCH --partition=knlall
#SBATCH --constraint knl,quad,cache
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=32
#SBATCH --output=p0.%j.%N.out
#SBATCH --error=p0.%j.%N.error
#SBATCH --mail-user=mraj@anl.gov
#SBATCH --time=01:00:00


# executable
exe=./ptrace


# max_step
dtSim=$((5000*7200*10))

# stepsize
dtParticle=$((300000/6))


# with prediction
prediction=0

# seed rate
seed_rate=1


# sample percentage
pred_percent=10



# num_procs=128
# seed_rate=16
# particle_file="particles/particles128.nc"
# gp_file="graphinfo/graph.info.part.128"
# opts="--blocks 128 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent"
# mpiexec -n $num_procs $exe $args

# num_procs=256
# seed_rate=8
# particle_file="particles/particles256.nc"
# gp_file="graphinfo/graph.info.part.256"
# opts="--blocks 256 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent"
# mpiexec -n $num_procs $exe $args

# num_procs=512
# seed_rate=4
# particle_file="particles/particles512.nc"
# gp_file="graphinfo/graph.info.part.512"
# opts="--blocks 512 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent"
# mpiexec -n $num_procs $exe $args

# num_procs=1024
# seed_rate=2
# particle_file="particles/particles1024.nc"
# gp_file="graphinfo/graph.info.part.1024"
# opts="--blocks 1024 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent"
# mpiexec -n $num_procs $exe $args


num_procs=2048
seed_rate=1
particle_file="particles/particles2048.nc"
gp_file="graphinfo/graph.info.part.2048"
opts="--blocks 2048 --check 0"
args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent"
mpiexec -n $num_procs $exe $args






