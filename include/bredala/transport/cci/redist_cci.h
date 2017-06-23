//---------------------------------------------------------------------------
//
// data interface
//
// Matthieu Dreher
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// mdreher@anl.gov
//
//--------------------------------------------------------------------------

#ifndef DECAF_REDIST_CCI_HPP
#define DECAF_REDIST_CCI_HPP

#include <iostream>
#include <assert.h>
#include <string>
#include <map>
#include <vector>

#include <bredala/transport/redist_comp.h>
#include <bredala/transport/mpi/types.h>
#include <cci.h>
#include <bredala/transport/cci/msg.h>

namespace decaf
{


class RedistCCI : public RedistComp
{
public:

    RedistCCI() :
        task_communicator_(MPI_COMM_NULL){}
    RedistCCI(int rankSource,
              int nbSources,
              int rankDest,
              int nbDests,
              int global_rank,
              CommHandle communicator,
              std::string name,
              RedistCommMethod commMethod,
              MergeMethod mergeMethod);
    virtual ~RedistCCI();

    virtual void flush();
    void shutdown();

    virtual void clearBuffers();

protected:

    // Compute the values necessary to determine how the data should be
    // splitted and redistributed.
    virtual void computeGlobal(pConstructData& data, RedistRole role) = 0;

    // Seperate the Data into chunks for each destination involve in the
    // component and fill the splitChunks vector
    virtual void splitData(pConstructData& data, RedistRole role) = 0;

    // Seperate system only data model. The data won't be split but duplicated
    void splitSystemData(pConstructData& data, RedistRole role);

    // Transfer the chunks from the sources to the destination. The data should be
    // be stored in the vector receivedChunks
    virtual void redistribute(pConstructData& data, RedistRole role);

    //Redistribution method using only point to point All to All communication
    void redistributeP2P(pConstructData& data, RedistRole role);

    //Redistribution method collective and sending only the necessary messages
    void redistributeCollective(pConstructData& data, RedistRole role);

    void init_connection_client(std::vector<std::string>&  server_uris);
    void init_connection_server(int nb_connections);

    CommHandle task_communicator_;      // communicator for all the processes involved in redist
    int task_rank_;                     // Rank of the process within the local task
    int task_size_;                     // Number of processes within the local task
    //std::vector<CommRequest> reqs;    // pending communication requests
    pConstructData transit;             // used when a source and destination are overlapping
    int *sum_;                          // used by the producer
    int *destBuffer_;
    std::string name_;                                  // Used to identify the file to exchange the URIs of the consumer
    cci_endpoint_t* endpoint_con_;                      // Endpoint of a consumer
    cci_endpoint_t* endpoint_prod_;                     // Endpoint of a producer
    std::map<cci_connection_t *, CCIchannel> channels_con_;     // Connections of the consumer side
    std::vector<CCIchannel> channels_prod_;
    std::map<cci_connection_t *, unsigned int> indexes_prod_;
    std::vector<cci_event_t *> event_queue_con_;                // Stack of events to process.
    std::vector<cci_event_t *> event_queue_prod_;

    std::vector< pConstructData > splitBuffer_;	// Buffer of container to avoid reallocation// used by the consumer
};

} // namespace

#endif
