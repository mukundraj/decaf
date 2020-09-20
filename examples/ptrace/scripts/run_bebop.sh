#!/bin/bash

#SBATCH --job-name=p1
#SBATCH --account=pedal
#SBATCH --partition=knlall
#SBATCH --constraint knl,quad,cache
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --output=p0.%j.%N.out
#SBATCH --error=p0.%j.%N.error
#SBATCH --mail-user=mraj@anl.gov
#SBATCH --time=01:00:00


# executable
exe=./ptrace


# max_step
# dtSim=$((5000*7200*10))
dtSim=$((5000))

# stepsize
dtParticle=$((300000/6))

# seed rate
seed_rate=1


# sample percentage
pred_percent=10


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

# num_procs=16
# particle_file="particles/particles16.nc"
# gp_file="graphinfo/graph.info.part.16"
# opts="--blocks 16 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle"
# mpiexec -n $num_procs $exe $args

num_procs=32
particle_file="particles/particles32.nc"
gp_file="graphinfo/graph.info.part.32"
opts="--blocks 32 --check 0"
args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle"
mpiexec -n $num_procs $exe $args

# num_procs=64
# particle_file="particles/particles64.nc"
# gp_file="graphinfo/graph.info.part.64"
# opts="--blocks 64 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle"
# mpiexec -n $num_procs $exe $args

# num_procs=128
# particle_file="particles/particles128.nc"
# gp_file="graphinfo/graph.info.part.128"
# opts="--blocks 128 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent"
# mpiexec -n $num_procs $exe $args

# num_procs=256
# particle_file="particles/particles256.nc"
# gp_file="graphinfo/graph.info.part.256"
# opts="--blocks 256 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent"
# mpiexec -n $num_procs $exe $args

# num_procs=512
# particle_file="particles/particles512.nc"
# gp_file="graphinfo/graph.info.part.512"
# opts="--blocks 512 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent"
# mpiexec -n $num_procs $exe $args

# num_procs=1024
# particle_file="particles/particles1024.nc"
# gp_file="graphinfo/graph.info.part.1024"
# opts="--blocks 1024 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent"
# mpiexec -n $num_procs $exe $args

# num_procs=2048
# particle_file="particles/particles2048.nc"
# gp_file="graphinfo/graph.info.part.2048"
# opts="--blocks 2048 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent"
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


# num_procs=16
# particle_file="particles/particles16.nc"
# gp_file="graphinfo/graph.info.part.16"
# opts="--blocks 16 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle"
# mpiexec -n $num_procs $exe $args

# num_procs=32
# particle_file="particles/particles32.nc"
# gp_file="graphinfo/graph.info.part.32"
# opts="--blocks 32 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle"
# mpiexec -n $num_procs $exe $args

# num_procs=64
# particle_file="particles/particles64.nc"
# gp_file="graphinfo/graph.info.part.64"
# opts="--blocks 64 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent"
# mpiexec -n $num_procs $exe $args

# num_procs=128
# particle_file="particles/particles128.nc"
# gp_file="graphinfo/graph.info.part.128"
# opts="--blocks 128 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent"
# mpiexec -n $num_procs $exe $args

# num_procs=256
# particle_file="particles/particles256.nc"
# gp_file="graphinfo/graph.info.part.256"
# opts="--blocks 256 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent"
# mpiexec -n $num_procs $exe $args

# num_procs=512
# particle_file="particles/particles512.nc"
# gp_file="graphinfo/graph.info.part.512"
# opts="--blocks 512 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent"
# mpiexec -n $num_procs $exe $args

# num_procs=1024
# particle_file="particles/particles1024.nc"
# gp_file="graphinfo/graph.info.part.1024"
# opts="--blocks 1024 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent"
# mpiexec -n $num_procs $exe $args

# num_procs=2048
# particle_file="particles/particles2048.nc"
# gp_file="graphinfo/graph.info.part.2048"
# opts="--blocks 2048 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent"
# mpiexec -n $num_procs $exe $args