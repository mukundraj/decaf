# converts an nx graph into a workflow data structure and runs the workflow

import imp
import sys
import os
import exceptions
import getopt
import argparse

import json
from collections import defaultdict

import fractions  #Definition of the least common multiple of two numbers
def lcm(a,b): return abs(a * b) / fractions.gcd(a,b) if a and b else 0

""" Object holding information about input/output data and types """
class Contract:
  """ Object holder for informations about the plateform to use """

  def __init__(self, jsonfilename = ""):
    self.inputs = {}
    self.outputs = {}
    if(jsonfilename != ""):
      with open(jsonfilename, "r") as f:
            content = f.read()
      contract = json.loads(content)
      if "inputs" in contract:
        for key, val in contract["inputs"].items():
          if len(val) == 1:
            val.append(1)
          self.inputs[key] = val
      if "outputs" in contract:
        for key, val in contract["outputs"].items():
          if len(val) == 1:
            val.append(1)
          self.outputs[key] = val
  # End constructor

  def addInput(self, key, typename, period=1):
    self.inputs[key] = [typename, period]

  def addInputFromDict(self, dict):
    for key, val in dict.items():
      if len(val) == 1:
        val.append(1) # period=1 by default
      self.inputs[key] = val

  def addOutput(self, key, typename, period=1):
    self.outputs[key] = [typename, period]

  def addOutputFromDict(self, dict):
    for key, val in dict.items():
      if len(val) == 1:
        val.append(1) # period=1 by default
      self.outputs[key] = val

# End of class Contract

# Checks if the intersection of each pair of prod/con contracts of an edge is non empty
# Then checks if all keys of a consumer are received
def check_contracts(graph, check_types):
    dict = defaultdict(set) # dictionary that lists the keys received by a consumer
    
    for edge in graph.edges_iter(data=True):
        prod = graph.node[edge[0]]
        con = graph.node[edge[1]]

        if ("contract" in prod) and ("contract" in con):
          if prod["contract"].outputs == {}:
              raise ValueError("ERROR while checking contracts: %s has no outputs" % edge[0])
          else:
              prod_out = prod["contract"].outputs

          if con["contract"].inputs == {}:
              raise ValueError("ERROR while checking contracts: %s has no inputs" % edge[1])
          else:
              con_in = con["contract"].inputs

          # We add for each edge the list of pairs key/type which is the intersection of the two contracts
          intersection_keys = {key:[con_in[key][0]] for key in con_in.keys() if (key in prod_out) and ( (check_types == 0) or (con_in[key][0] == prod_out[key][0]) ) }
          for key in intersection_keys.keys():
              dict[edge[1]].add(key)
              intersection_keys[key].append(lcm(con_in[key][1], prod_out[key][1]))

          if len(intersection_keys) == 0:
              raise ValueError("ERROR intersection of keys from %s and %s is empty" % (edge[0], edge[1]))

          # We add the list of pairs key/type to be exchanged b/w the producer and consumer
          edge[2]['keys'] = intersection_keys


    for name, set_k in dict.items():
    	needed = graph.node[name]['contract'].inputs
    	received = list(set_k)
        complement = [item for item in needed if item not in received]
        if len(complement) != 0:
          s = "ERROR %s does not receive the following:" % name
          for key in complement:
            s+= " %s:%s," %(key, needed[key])
          #print s.rstrip(',')
          raise ValueError(s.rstrip(','))
# End of function check_contracts


class Topology:
  """ Object holder for informations about the plateform to use """

  def __init__(self, name, nProcs = 0, hostfilename = "", offsetRank = 0, procPerNode = 0, offsetProcPerNode = 0):

      self.name = name                  # Name of the topology
      self.nProcs = nProcs                # Total number of MPI ranks in the topology
      self.hostfilename = hostfilename  # File containing the list of hosts. Must match  totalProc
      self.hostlist = []                # List of hosts.
      self.offsetRank = offsetRank      # Rank of the first rank in this topology
      self.nNodes = 0                   # Number of nodes
      self.nodes = []                   # List of nodes
      self.procPerNode = procPerNode    # Number of MPI ranks per node. If not given, deduced from the host list
      self.offsetProcPerNode = offsetProcPerNode    # Offset of the proc ID to start per node with process pinning


      if hostfilename != "":
        f = open(hostfilename, "r")
        content = f.read()
        content = content.rstrip(' \t\n\r') #Clean the last character, usually a \n
        self.hostlist = content.split('\n')
        self.nodes = list(set(self.hostlist))   # Removing the duplicate hosts
        self.nNodes = len(self.nodes)
      else:
        self.hostlist = ["localhost"] * self.nProcs # should be the same as the 3 above lines
        self.nodes = ["localhost"]
        self.nNodes = 1

      if len(self.hostlist) != self.nProcs:
          raise ValueError("The number of hosts does not match the number of procs.")
  # End of constructor

  def isInitialized(self):
      return self.nProcs > 0


  def toStr(self):
      content = ""
      content += "name: " + self.name + "\n"
      content += "nProcs: " + str(self.nProcs) + "\n"
      content += "offsetRank: " + str(self.offsetRank) + "\n"
      content += "host list: " + str(self.hostlist) + "\n"

      return content

  def subTopology(self, name, nProcs, procOffset):

      # Enough rank check
      if nProcs > self.nProcs:
        raise ValueError("Not enough rank available")

      subTopo = Topology(name, nProcs, offsetRank = procOffset)

      # Adjusting the hostlist
      subTopo.hostlist = self.hostlist[procOffset:procOffset+nProcs]
      subTopo.nodes = list(set(subTopo.hostlist))
      subTopo.nNodes = len(subTopo.nodes)

      return subTopo

  def splitTopology(self, names, nProcs):

      # Same size list check
      if len(nProcs) != len (names):
        raise ValueError("The size of the list of procs and names don't match.")

      # Enough rank check
      if sum(nProcs) > self.nProcs:
        raise ValueError("Not enough rank available")

      offset = 0
      splits = []

      for i in range(0, len(names)):
        subTopo = Topology(names[i], nProcs[i], offsetRank = offset)

        subTopo.hostlist = self.hostlist[offset:offset+nProcs[i]]
        subTopo.nodes = list(set(subTopo.hostlist))
        subTopo.nNodes = len(subTopo.nodes)
        offset += nProcs[i]

        splits.append(subTopo)

      return splits
# End of class Topology


def initParserForTopology(parser):
    """ Add the necessary arguments for initialize a topology.
        The parser might be completed by the user in other functions
    """

    parser.add_argument(
        "-n", "--np",
        help = "Total number of MPI ranks used in the workflow",
        type=int,
        default=1,
        required=True,
        dest="nprocs"
        )
    parser.add_argument(
        "--hostfile",
        help = "List of host for each MPI rank. Must match the number of MPI ranks.",
        default="",
        dest="hostfile"
        )

def topologyFromArgs(args):
    return Topology("global", args.nprocs, hostfilename = args.hostfile)

def processTopology(graph):
    """ Check all nodes and edge if a topology is present.
        If yes, fill the fields start_proc and nprocs """

    for node in graph.nodes_iter(data=True):
        if 'topology' in node[1]:
            topo = node[1]['topology']
            node[1]['start_proc'] = topo.offsetRank
            node[1]['nprocs'] = topo.nProcs

    for edge in graph.edges_iter(data=True):
        if 'topology' in edge[2]:
            topo = edge[2]['topology']
            edge[2]['start_proc'] = topo.offsetRank
            edge[2]['nprocs'] = topo.nProcs


def workflowToJson(graph, outputFile, check_types):
    print "Generating graph description file "+outputFile

    nodes   = []
    links   = []

    f = open(outputFile, "w")
    content = ""
    content +="{\n"
    content +="   \"workflow\" :\n"
    content +="   {\n"
    content +="   \"check_types\" : \""+str(check_types)+"\",\n"
    content +="   \"nodes\" : [\n"

    i = 0
    for node in graph.nodes_iter(data=True):
        content +="      {\n"
        content +="       \"start_proc\" : "+str(node[1]['start_proc'])
        content +=",\n       \"nprocs\" : "+str(node[1]['nprocs'])
        content +=",\n       \"func\" : \""+node[1]['func']+"\""
        content +="\n      },\n"

        node[1]['index'] = i
        i += 1

    content = content.rstrip(",\n")
    content += "\n"
    content +="   ],\n"
    content +="\"edges\" : [\n"

    for edge in graph.edges_iter(data=True):
        prod  = graph.node[edge[0]]['index']
        con   = graph.node[edge[1]]['index']
        content +="      {\n"
        content +="       \"start_proc\" : "+str(edge[2]['start_proc'])
        content +=",\n       \"nprocs\" : "+str(edge[2]['nprocs'])
        content +=",\n       \"source\" : "+str(prod)
        content +=",\n       \"target\" : "+str(con)
        content +=",\n       \"prod_dflow_redist\" : \""+edge[2]['prod_dflow_redist']+"\""

        if(edge[2]['nprocs'] != 0):
          content +=",\n       \"dflow_con_redist\" : \""+edge[2]['dflow_con_redist']+"\""
          content +=",\n       \"func\" : \""+edge[2]['func']+"\""
          content +=",\n       \"path\" : \""+edge[2]['path']+"\""

        if "keys" in edge[2]:
          content +=",\n       \"keys\" : "+json.dumps(edge[2]['keys'], sort_keys=True)

        content +="\n      },\n"

    content = content.rstrip(",\n")
    content += "\n"
    content +="   ]\n"
    content +="   }\n"
    content +="}"

    f.write(content)
    f.close()

# Looking for a node/edge starting at a particular rank
def getNodeWithRank(rank, graph):

    for node in graph.nodes_iter(data=True):
        if node[1]['start_proc'] == rank:
            return ('node', node)

    for edge in graph.edges_iter(data=True):
        if (edge[2]['start_proc'] == rank) and (edge[2]['nprocs'] != 0):
            return ('edge', edge)

    return ('notfound', graph.nodes_iter())


# Build the mpirun command in MPMD mode
# Parse the graph to sequence the executables with their
# associated MPI ranks and arguments
def workflowToSh(graph, outputFile, mpirunOpt = "", mpirunPath = ""):
    print "Generating bash command script "+outputFile

    currentRank = 0
    mpirunCommand = mpirunPath+"mpirun "+mpirunOpt+" "
    nbExecutables = graph.number_of_nodes() + graph.number_of_edges()

    # Parsing the graph looking at the current rank
    for i in range(0, nbExecutables):
      (type, exe) = getNodeWithRank(currentRank, graph)
      if type == 'none':
        print 'ERROR: Unable to find an executable for the rank ' + str(rank)
        exit()

      if type == 'node':

        if exe[1]['nprocs'] == 0:
          print 'ERROR: a node can not have 0 MPI rank.'
          exit()

        mpirunCommand += "-np "+str(exe[1]['nprocs'])
        mpirunCommand += " "+str(exe[1]['cmdline'])+" : "
        currentRank += exe[1]['nprocs']

      if type == 'edge':

        if exe[2]['nprocs'] != 0:
          mpirunCommand += "-np "+str(exe[2]['nprocs'])
          mpirunCommand += " "+str(exe[2]['cmdline'])+" : "
          currentRank += exe[2]['nprocs']

    # Writting the final file
    f = open(outputFile, "w")
    content = ""
    content +="#! /bin/bash\n\n"
    content +=mpirunCommand
    content = content.rstrip(" : ")

    f.write(content)
    f.close()
    os.system("chmod a+rx "+outputFile)

""" Process the graph and generate the necessary files
  check_types = 0 means no typechecking; 1 means typechecking at python script only
  2 means typechecking both at python script and at runtime; Only relevant if contracts are used
"""
def processGraph(graph, name, check_types = 1, mpirunPath = "", mpirunOpt = ""):
    check_contracts(graph, check_types)
    processTopology(graph)
    workflowToJson(graph, name+".json", check_types)
    workflowToSh(graph, name+".sh", mpirunOpt = mpirunOpt, mpirunPath = mpirunPath)
