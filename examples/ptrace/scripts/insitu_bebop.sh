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
#SBATCH --time=00:15:00

export DECAF_FOLDER=/home/mraj/tests/mpas-decaf-ptrace/decaf/
export BOOST_INSTALL_DIR=/home/mraj/libraries/boost_1_65_1/install/

export NETCDF=/home/mraj/libs_mpas/install
export PNETCDF=/home/mraj/libs_mpas/install
export PIO=/home/mraj/libs_mpas/install

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DECAF_FOLDER/install/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BOOST_INSTALL_DIR/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/mraj/libs_mpas/mpich-3.2.1/install/lib/

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

# skip_rate=1
# particle_file="/home/mraj/tests/parfiles/particles8.nc"
# gp_file="graphinfo/graph.info.part.8"
# opts="--blocks 8 --check 0"
# args="$opts $particle_file $gp_file $prediction $dtSim $dtParticle $seed_rate $pred_percent $skip_rate"

# mpiexec -n $num_procs $exe $args



#srun -n 16 -l --multi-prog my.config08
srun -n 64 -l  --multi-prog my.config32
