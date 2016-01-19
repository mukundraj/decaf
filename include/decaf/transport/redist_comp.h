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

#ifndef DECAF_REDIST_COMP_H
#define DECAF_REDIST_COMP_H

#include <decaf/data_model/basedata.h>
#include <vector>
#include <memory>
#include <iostream>

namespace decaf
{

    enum RedistRole
    {
        DECAF_REDIST_SOURCE,
        DECAF_REDIST_DEST,
    };

    enum mpiTags {
        MPI_METADATA_TAG = 1,
        MPI_DATA_TAG,
    };

    // This class defines the common interface for the redistribution component (MxN)
    // This interface is independant from the datatype or the transport
    // implementation. Specialized components will implement the redistribution
    // in fonction of transport layer
    class RedistComp
    {
    public:
        RedistComp(){}
        RedistComp(int rankSource,
                   int nbSources,
                   int rankDest,
                   int nbDests) :
            rankSource_(rankSource),
            nbSources_(nbSources),
            rankDest_(rankDest),
            nbDests_(nbDests),
            summerizeDest_(NULL) {}

        virtual ~RedistComp(){}

        // Run the pipeline of operations to redistribute the data.
        // This function is the only one that should be called from
        // the main programm

        /* virtual int process(std::shared_ptr<BaseData> data, RedistRole role); */
        int process(std::shared_ptr<BaseData> data, RedistRole role);

        int getRankSource();
        int getNbSources();
        int getRankDest();
        int getNbDest();

        bool isSource();
        bool isDest();

        virtual void flush() = 0;

    protected:

        // Compute the values necessary to determine how the data should be
        // splitted and redistributed.
        virtual void computeGlobal(std::shared_ptr<BaseData> data, RedistRole role)=0;

        // Seperate the Data into chunks for each destination involve in the
        // component and fill the splitChunks vector
        virtual void splitData(std::shared_ptr<BaseData> data, RedistRole role)=0;

        // Transfer the chunks from the sources to the destination. The data
        // should be stored in the vector receivedChunks
        virtual int redistribute(std::shared_ptr<BaseData> data, RedistRole role)=0;

        // Merge the chunks from the vector receivedChunks into one single Data.
        std::shared_ptr<BaseData> merge(RedistRole role);

        int rankSource_; // Rank of the first source (=sender)
        int nbSources_;  // Number of sources, supposed to be consecutives
        int rankDest_;   // Rank of the first destination (=receiver)
        int nbDests_;     // Number of destinations

        // TODO: check for redundancy from above members
        int rank_;                // Rank in the group communicator
        int size_;                // Size of the group communicator
        int local_source_rank_;   // Rank of the first source in communicator_
        int local_dest_rank_;     // Rank of the first destination in communicator_

        std::vector<std::shared_ptr<BaseData> > splitChunks_;
        std::vector<std::shared_ptr<char> > receivedChunks_;
        int* summerizeDest_;
        std::vector<int> destList_;
    };

} //namespace decaf

#endif
