//---------------------------------------------------------------------------
//
// workflow definition
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
//--------------------------------------------------------------------------

#ifndef DECAF_WORKFLOW_HPP
#define DECAF_WORKFLOW_HPP

#include "decaf.hpp"

#include <stdio.h>
#include <vector>
#include <queue>
#include <string>
#include <iostream>
#include <fstream>

#include "json/json.h"

using namespace std;

struct WorkflowNode                          // a producer or consumer
{
    WorkflowNode()                                {}
    WorkflowNode(int start_proc_,
                 int nprocs_,
                 string func_) :
        start_proc(start_proc_),
        nprocs(nprocs_),
        func(func_),
        args(NULL)                                {}
    vector<int> out_links; // indices of outgoing links
    vector<int> in_links;  // indices of incoming links
    int start_proc;        // starting proc rank in world communicator for this producer or consumer
    int nprocs;            // number of processes for this producer or consumer
    string func;           // name of node callback
    void* args;            // callback arguments
    void add_out_link(int link) { out_links.push_back(link); }
    void add_in_link(int link) { in_links.push_back(link); }
};

struct WorkflowLink                          // a dataflow
{
    WorkflowLink()                                {}
    WorkflowLink(int prod_,
                 int con_,
                 int start_proc_,
                 int nprocs_,
                 string func_,
                 string path_,
                 string prod_dflow_redist_,
                 string dflow_con_redist_) :
        prod(prod_),
        con(con_),
        start_proc(start_proc_),
        nprocs(nprocs_),
        func(func_),
        args(NULL),
        path(path_),
        prod_dflow_redist(prod_dflow_redist_),
        dflow_con_redist(dflow_con_redist_)       {}
    int prod;                   // index in vector of all workflow nodes of producer
    int con;                    // index in vector of all workflow nodes of consumer
    int start_proc;             // starting process rank in world communicator for the dataflow
    int nprocs;                 // number of processes in the dataflow
    string func;                // name of dataflow callback
    void* args;                 // callback arguments
    string path;                // path to callback function module
    string prod_dflow_redist;   // redistribution component between producer and dflow
    string dflow_con_redist;    // redistribution component between dflow and consumer
};

struct Workflow                              // an entire workflow
{
    Workflow()                                    {}
    Workflow(vector<WorkflowNode>& nodes_,
             vector<WorkflowLink>& links_) :
        nodes(nodes_),
        links(links_)                             {}
    vector<WorkflowNode> nodes;              // all the workflow nodes
    vector<WorkflowLink> links;              // all the workflow links
    bool my_node(int proc, int node)         // whether my process is part of this node
        {
            return(proc >= nodes[node].start_proc &&
                   proc <  nodes[node].start_proc + nodes[node].nprocs);
        }
    bool my_link(int proc, int link)         // whether my process is part of this link
        {
            return(proc >=links[link].start_proc &&
                   proc < links[link].start_proc + links[link].nprocs);
        }
    bool my_in_link(int proc, int link)      // whether my process gets input data from this link
        {
            for (size_t i = 0; i< nodes.size(); i++)
            {
                if (proc >= nodes[i].start_proc && // proc is mine
                    proc <  nodes[i].start_proc + nodes[i].nprocs)
                {
                    for (size_t j = 0; j < nodes[i].in_links.size(); j++)
                    {
                        if (nodes[i].in_links[j] == link)
                            return true;
                    }
                }
            }
            return false;
        }
    bool my_out_link(int proc, int link)      // whether my process puts output data to this link
        {
            for (size_t i = 0; i< nodes.size(); i++)
            {
                if (proc >= nodes[i].start_proc && // proc is mine
                    proc <  nodes[i].start_proc + nodes[i].nprocs)
                {
                    for (size_t j = 0; j < nodes[i].out_links.size(); j++)
                    {
                        if (nodes[i].out_links[j] == link)
                            return true;
                    }
                }
            }
            return false;
        }

  static void
  make_wflow_from_json( Workflow& workflow, string json_path )
  {
    char *prefix = getenv("DECAF_PREFIX");
    if( prefix == nullptr ) {
      fprintf( stderr, "ERROR: environment variable DECAF_PREFIX not defined."
	       "Please export DECAF_PREFIX to point to the root of your decaf"
	       "install directory.\n");
      exit(1);
    }

    Json::Value root;

    try {
      std::ifstream fs( json_path, std::ifstream::binary );
      fs >> root;
    }
    catch( std::exception& e ) {
      cerr << "std::exception while opening/parsing JSON file("
	   << json_path
	   << ")"
	   << e.what()
	   << endl;
    }

    Json::Value nodes = root["nodes"];
    Json::Value edges = root["edges"];

    for( auto &&v : nodes ) {
      /* 
       * iterate over the list of nodes, creating and populating WorkflowNodes as we go
       */
      WorkflowNode node;
      node.out_links.clear();
      node.in_links.clear();
      /* we defer actually linking nodes until we read the edge list */
      
      node.start_proc = v.get("start_proc").asInt();
      node.nprocs = v.get("nprocs").asInt();
      node.func = v.get("func").asString();
      workflow.nodes.push_back( node );
    }

    string path = string( prefix, strlen(prefix) );
    path.append( root["path"].asString() );

    for( auto &&v : edges ) {

      WorkflowLink link;
      
      /* link the nodes correctly(?) */      
      link.prod = v.get("source").asInt();
      link.con = v.get("target").asInt();
      workflow.nodes.at( link.prod ).out_links.push_back( link.con );
      workflow.nodes.at( link.con ).in_links.push_back( link.prod );

      link.path = path;
      link.start_proc = v.get("start_proc").asInt();
      link.nprocs = v.get("nprocs").asInt();
      link.func = v.get("func").asString();
      link.prod_dflow_redist = v.get("count").asString();
      link.dflow_con_redist = v.get("count").asString();
      workflow.links.push_back( link );
    }
  }
  
};

#endif
