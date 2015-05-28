# a small 2-node example, just a producer and consumer

import networkx as nx
import os

# --- set your options here ---

# path to .so driver
driver_path = os.environ['DECAF_PREFIX'] + '/examples/direct/python/libpy_direct.so'

# path to .so module for callback functions
mod_path = os.environ['DECAF_PREFIX'] + '/examples/direct/libmod_direct.so'

# define workflow graph
# 2-node workflow
#
#    lammps (4 procs) - print (2 procs)
#
#  entire workflow takes 4 procs because of the overlap

w = nx.DiGraph()
w.add_node("prod", start_proc=0, nprocs=4, prod_func='prod', con_func=''       , path=mod_path)
w.add_node("con",  start_proc=2, nprocs=2, prod_func= ''   , con_func='con'    , path=mod_path)
w.add_edge("prod", "con", start_proc=1, nprocs=2,            dflow_func='dflow', path=mod_path)

# total number of time steps
prod_nsteps  = 2
con_nsteps   = 2

# --- do not edit below this point --

import imp
driver = imp.load_dynamic('driver', driver_path)
driver.pyrun(w, prod_nsteps, con_nsteps)
