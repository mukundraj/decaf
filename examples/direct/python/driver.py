# a small 2-node example, just a producer and consumer

import networkx as nx

# --- set your options here ---

# path to .so module
#path = '/Users/tpeterka/software/decaf/install/examples/direct/python/libpy_direct.so'
path = '/home/matthieu/Argonne/Decaf/build/examples/direct/python/libpy_direct.so'

# define workflow graph
# 2-node workflow
#
#    lammps (4 procs) - print (2 procs)
#
#  entire workflow takes 8 procs (2 dataflow procs between producer and consumer)
#  dataflow can be overlapped, but currently all disjoint procs (simplest case)

w = nx.DiGraph()

# example of 4 nodes and 3 edges (single source)
# this is the example diagrammed above, and for which driver.pyx is made
w.add_node("prod", start_proc=0, nprocs=4, prod_func='prod'     , con_func=''        )
w.add_node("con",  start_proc=6, nprocs=2, prod_func= ''        , con_func='con'     )
w.add_edge("prod", "con", start_proc=4, nprocs=2, dflow_func='dflow'                 )

# total number of time steps
prod_nsteps  = 1
con_nsteps   = 1

# --- do not edit below this point --

import imp
driver = imp.load_dynamic('driver', path)
driver.pyrun(w, prod_nsteps, con_nsteps)
