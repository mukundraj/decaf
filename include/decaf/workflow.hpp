﻿//---------------------------------------------------------------------------
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

//#include "decaf.hpp"

#include <stdio.h>
#include <vector>
#include <queue>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"

#include <manala/types.h>
#include <decaf/tools.hpp>
#include <manala/tools.h>

namespace bpt = boost::property_tree;

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
    void add_out_link(int link);
    void add_in_link(int link);
};

struct WorkflowLink                          // a dataflow
{
    WorkflowLink()                                {}
    /*WorkflowLink(int prod_,			// This constructor is never used
                                 int con_,
                 int start_proc_,
                 int nprocs_,
                 string func_,
                 string path_,
                 string prod_dflow_redist_,
                                 string dflow_con_redist_,
                                 vector<ContractKey> list_keys_,
                                 Check_level check_level_,
                 string stream_) :
        prod(prod_),
        con(con_),
        start_proc(start_proc_),
        nprocs(nprocs_),
        func(func_),
        args(NULL),
        path(path_),
        prod_dflow_redist(prod_dflow_redist_),
                dflow_con_redist(dflow_con_redist_),
                list_keys(list_keys_),
                check_level(check_level_),
        stream(stream_){} */
    int prod;                   // index in vector of all workflow nodes of producer
    int con;                    // index in vector of all workflow nodes of consumer
    int start_proc;             // starting process rank in world communicator for the dataflow
    int nprocs;                 // number of processes in the dataflow
    string func;                // name of dataflow callback
    void* args;                 // callback arguments
    string path;                // path to callback function module
    string prod_dflow_redist;   // redistribution component between producer and dflow
    string dflow_con_redist;    // redistribution component between dflow and consumer
    /*string stream;              // Type of stream policy to use (none, single, double)
    string frame_policy;        // Policy to use to manage the incoming frames
    unsigned int prod_freq_output;              // Output frequency of the producer
    string storage_policy;                      // Type of storage collection to use
    vector<StorageType> storages;               // Different level of storage availables
    vector<unsigned int> storage_max_buffer;    // Maximum number of frame*/
    ManalaInfo manala_info;

    string srcPort;		// Portname of the source
    string destPort;		// Portname of the dest

    vector<ContractKey> keys_link;   // List of keys to be exchanged b/w the link and the consumer
    vector<ContractKey> list_keys;   // list of key to be exchanged b/w the producer and consumer or producer and link
    Check_level check_level;		 // level of checking for the types of data to be exchanged
    bool bAny;						 // Whether the filtering will check the contracts but keep any other field or not


};

struct Workflow                              // an entire workflow
{
    Workflow()                                    {}
    Workflow(vector<WorkflowNode>& nodes_,
             vector<WorkflowLink>& links_) :
        nodes(nodes_),
        links(links_)                             {}
    vector<WorkflowNode> nodes;             // all the workflow nodes
    vector<WorkflowLink> links;             // all the workflow links
    bool my_node(int proc, int node);       // whether my process is part of this node

    bool my_link(int proc, int link);       // whether my process is part of this link

    bool my_in_link(int proc, int link);    // whether my process gets input data from this link

    bool my_out_link(int proc, int link);   // whether my process puts output data to this link


    static void
    make_wflow_from_json( Workflow& workflow, const string& json_path );
};

#endif
