#!/bin/bash

#SBATCH --job-name=p1
#SBATCH --account=pedal
#SBATCH --partition=knlall
#SBATCH --constraint knl,quad,cache
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --output=p0.%j.%N.out
#SBATCH --error=p0.%j.%N.error
#SBATCH --mail-user=mraj@anl.gov
#SBATCH --time=01:00:00


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
particle_file="particles/particles32.nc"
gp_file="graphinfo/graph.info.part.32"
opts="--blocks 32 --check 0"
args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent $skip_rate"
mpiexec -n $num_procs $exe $args


# skip_rate=2
# particle_file="particles/particles512.nc"
# gp_file="graphinfo/graph.info.part.512"
# opts="--blocks 512 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent $skip_rate"
# mpiexec -n $num_procs $exe $args

# skip_rate=5
# particle_file="particles/particles512.nc"
# gp_file="graphinfo/graph.info.part.512"
# opts="--blocks 512 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent $skip_rate"
# mpiexec -n $num_procs $exe $args


# skip_rate=10
# particle_file="particles/particles512.nc"
# gp_file="graphinfo/graph.info.part.512"
# opts="--blocks 512 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent $skip_rate"
# mpiexec -n $num_procs $exe $args

# skip_rate=20
# particle_file="particles/particles512.nc"
# gp_file="graphinfo/graph.info.part.512"
# opts="--blocks 512 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent $skip_rate"
# mpiexec -n $num_procs $exe $args

# skip_rate=50
# particle_file="particles/particles512.nc"
# gp_file="graphinfo/graph.info.part.512"
# opts="--blocks 512 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent $skip_rate"
# mpiexec -n $num_procs $exe $args

# skip_rate=100
# particle_file="particles/particles512.nc"
# gp_file="graphinfo/graph.info.part.512"
# opts="--blocks 512 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent $skip_rate"
# mpiexec -n $num_procs $exe $args
