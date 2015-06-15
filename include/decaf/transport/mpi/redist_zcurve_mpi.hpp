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

#ifndef DECAF_REDIST_ZCURVE_MPI_HPP
#define DECAF_REDIST_ZCURVE_MPI_HPP

#include <iostream>
#include <assert.h>

#include <decaf/redist_comp.h>


namespace decaf
{

unsigned int Morton_3D_Encode_10bit( unsigned int index1, unsigned int index2, unsigned int index3 )
{ // pack 3 10-bit indices into a 30-bit Morton code
  index1 &= 0x000003ff;
  index2 &= 0x000003ff;
  index3 &= 0x000003ff;
  index1 |= ( index1 << 16 );
  index2 |= ( index2 << 16 );
  index3 |= ( index3 << 16 );
  index1 &= 0x030000ff;
  index2 &= 0x030000ff;
  index3 &= 0x030000ff;
  index1 |= ( index1 << 8 );
  index2 |= ( index2 << 8 );
  index3 |= ( index3 << 8 );
  index1 &= 0x0300f00f;
  index2 &= 0x0300f00f;
  index3 &= 0x0300f00f;
  index1 |= ( index1 << 4 );
  index2 |= ( index2 << 4 );
  index3 |= ( index3 << 4 );
  index1 &= 0x030c30c3;
  index2 &= 0x030c30c3;
  index3 &= 0x030c30c3;
  index1 |= ( index1 << 2 );
  index2 |= ( index2 << 2 );
  index3 |= ( index3 << 2 );
  index1 &= 0x09249249;
  index2 &= 0x09249249;
  index3 &= 0x09249249;
  return( index1 | ( index2 << 1 ) | ( index3 << 2 ) );
}

void Morton_3D_Decode_10bit( const unsigned int morton,
                                 unsigned int& index1, unsigned int& index2, unsigned int& index3 )
    { // unpack 3 10-bit indices from a 30-bit Morton code
      unsigned int value1 = morton;
      unsigned int value2 = ( value1 >> 1 );
      unsigned int value3 = ( value1 >> 2 );
      value1 &= 0x09249249;
      value2 &= 0x09249249;
      value3 &= 0x09249249;
      value1 |= ( value1 >> 2 );
      value2 |= ( value2 >> 2 );
      value3 |= ( value3 >> 2 );
      value1 &= 0x030c30c3;
      value2 &= 0x030c30c3;
      value3 &= 0x030c30c3;
      value1 |= ( value1 >> 4 );
      value2 |= ( value2 >> 4 );
      value3 |= ( value3 >> 4 );
      value1 &= 0x0300f00f;
      value2 &= 0x0300f00f;
      value3 &= 0x0300f00f;
      value1 |= ( value1 >> 8 );
      value2 |= ( value2 >> 8 );
      value3 |= ( value3 >> 8 );
      value1 &= 0x030000ff;
      value2 &= 0x030000ff;
      value3 &= 0x030000ff;
      value1 |= ( value1 >> 16 );
      value2 |= ( value2 >> 16 );
      value3 |= ( value3 >> 16 );
      value1 &= 0x000003ff;
      value2 &= 0x000003ff;
      value3 &= 0x000003ff;
      index1 = value1;
      index2 = value2;
      index3 = value3;
    }


  class RedistZCurveMPI : public RedistComp
  {
  public:
      RedistZCurveMPI() :
          communicator_(MPI_COMM_NULL),
          commSources_(MPI_COMM_NULL),
          commDests_(MPI_COMM_NULL) {}
      RedistZCurveMPI(int rankSource, int nbSources,
                     int rankDest, int nbDests,
                     CommHandle communicator,
                     std::vector<float> bBox = std::vector<float>(),
                     std::vector<int> slices = std::vector<int>());
      ~RedistZCurveMPI();

  protected:

      // Compute the values necessary to determine how the data should be splitted
      // and redistributed.
      virtual void computeGlobal(std::shared_ptr<BaseData> data, RedistRole role);

      // Seperate the Data into chunks for each destination involve in the component
      // and fill the splitChunks vector
      virtual void splitData(std::shared_ptr<BaseData> data, RedistRole role);

      // Transfert the chunks from the sources to the destination. The data should be
      // be stored in the vector receivedChunks
      virtual void redistribute(std::shared_ptr<BaseData> data, RedistRole role);


      // Merge the chunks from the vector receivedChunks into one single Data.
      virtual std::shared_ptr<BaseData> merge(RedistRole role);

      // Merge the chunks from the vector receivedChunks into one single data->
      //virtual BaseData* Merge();

      //bool isSource(){ return rank_ <  nbSources_; }
      //bool isDest(){ return rank_ >= rankDest_ - rankSource_; }

      bool isSource(){ return rank_ >= local_source_rank_ && rank_ < local_source_rank_ + nbSources_; }
      bool isDest(){ return rank_ >= local_dest_rank_ && rank_ < local_dest_rank_ + nbDests_; }

      CommHandle communicator_; // Communicator for all the processes involve
      CommHandle commSources_;  // Communicator of the sources
      CommHandle commDests_;    // Communicator of the destinations
      std::vector<CommRequest> reqs;       // pending communication requests

      int rank_;                // Rank in the group communicator
      int size_;                // Size of the group communicator
      int local_source_rank_;   // Rank of the first source in communicator_
      int local_dest_rank_;     // Rank of the first destination in communicator_

      // We keep these values so we can reuse them between 2 iterations
      int global_item_rank_;    // Index of the first item in the global array
      int global_nb_items_;     // Number of items in the global array

      std::shared_ptr<BaseData> transit; // Used then a source and destination are overlapping

      bool bBBox_;
      std::vector<float> bBox_;
      std::vector<int> slices_;
      std::vector<float> slicesDelta_;
      int indexes_per_dest_;
      int rankOffset_;


  };

} // namespace

decaf::
RedistZCurveMPI::RedistZCurveMPI(int rankSource, int nbSources,
                               int rankDest, int nbDests, CommHandle world_comm,
                                 std::vector<float> bBox,
                                 std::vector<int> slices ) :
    RedistComp(rankSource, nbSources, rankDest, nbDests),
    communicator_(MPI_COMM_NULL),
    commSources_(MPI_COMM_NULL),
    commDests_(MPI_COMM_NULL),
    bBBox_(false),
    bBox_(bBox),
    slices_(slices)
{
    std::cout<<"Generation of the Redist component with ["<<rankSource_<<","
            <<nbSources_<<","<<rankDest_<<","<<nbDests_<<"]"<<std::endl;
    MPI_Group group, groupRedist, groupSource, groupDest;
    int range[3];
    MPI_Comm_group(world_comm, &group);
    int world_rank, world_size;

    MPI_Comm_rank(world_comm, &world_rank);
    MPI_Comm_size(world_comm, &world_size);

    local_source_rank_ = 0;                     //Rank of first source in communicator_
    local_dest_rank_ = rankDest_ - rankSource_; //Rank of first destination in communucator_
    std::cout<<"Local destination rank : "<<local_dest_rank_<<std::endl;

    //Generation of the group with all the sources and destination
    range[0] = rankSource;
    range[1] = max(rankSource + nbSources - 1, rankDest + nbDests - 1);
    range[2] = 1;

    MPI_Group_range_incl(group, 1, &range, &groupRedist);
    MPI_Comm_create_group(world_comm, groupRedist, 0, &communicator_);
    MPI_Comm_rank(communicator_, &rank_);
    MPI_Comm_size(communicator_, &size_);
    std::cout<<"Rank in the Redist component : "<<rank_<<std::endl;

    //Generation of the group with all the sources
    //if(world_rank >= rankSource_ && world_rank < rankSource_ + nbSources_)
    if(isSource())
    {
        std::cout<<"Generating the source communicator"<<std::endl;
        range[0] = rankSource;
        range[1] = rankSource + nbSources - 1;
        range[2] = 1;
        MPI_Group_range_incl(group, 1, &range, &groupSource);
        MPI_Comm_create_group(world_comm, groupSource, 0, &commSources_);
        MPI_Group_free(&groupSource);
        int source_rank;
        MPI_Comm_rank(commSources_, &source_rank);
        std::cout<<"Source Rank in the Redist component : "<<source_rank<<std::endl;
    }

    //Generation of the group with all the Destinations
    //if(world_rank >= rankDest_ && world_rank < rankDest_ + nbDests_)
    if(isDest())
    {
        std::cout<<"Generating the destination communicator"<<std::endl;
        range[0] = rankDest;
        range[1] = rankDest + nbDests - 1;
        range[2] = 1;
        MPI_Group_range_incl(group, 1, &range, &groupDest);
        MPI_Comm_create_group(world_comm, groupDest, 0, &commDests_);
        MPI_Group_free(&groupDest);
        int dest_rank;
        MPI_Comm_rank(commDests_, &dest_rank);
        std::cout<<"Dest Rank in the Redist component : "<<dest_rank<<std::endl;

    }

    MPI_Group_free(&group);
    MPI_Group_free(&groupRedist);

    if(slices_.size() != 3)
    {
        std::cout<<"No slices given, will use 8,8,8"<<std::endl;
        slices_ = { 8,8,8 };
    }
    else
    {
        for(unsigned int i = 0; i < 3; i++)
        {
            if(slices_.at(i) < 1)
            {
                std::cout<<" ERROR : slices can't be inferior to 1. Switching the value to 1."<<std::endl;
                slices_.at(i) = 1;
            }
        }
    }

    //Computing the index ranges per destination
    //We can use -1 because the constructor makes sure that all items of slices_ are > 0
    int maxIndex = Morton_3D_Encode_10bit(slices_[0]-1, slices_[1]-1, slices_[2]-1);
    indexes_per_dest_ = maxIndex /  nbDests_;
    rankOffset_ = maxIndex %  nbDests_;

    // Checking the bounding box and updating the slicesDelta if possible
    if(bBox_.size() == 6)
    {
        slicesDelta_.resize(3);
        slicesDelta_[0] = (bBox_[3] - bBox_[0]) / (float)(slices_[0]);
        slicesDelta_[1] = (bBox_[4] - bBox_[1]) / (float)(slices_[1]);
        slicesDelta_[2] = (bBox_[5] - bBox_[2]) / (float)(slices_[2]);
    }
}

decaf::
RedistZCurveMPI::~RedistZCurveMPI()
{
  if (communicator_ != MPI_COMM_NULL)
    MPI_Comm_free(&communicator_);
  if (commSources_ != MPI_COMM_NULL)
    MPI_Comm_free(&commSources_);
  if (commDests_ != MPI_COMM_NULL)
    MPI_Comm_free(&commDests_);
}

/*void
decaf::
RedistZCurveMPI::put(std::shared_ptr<BaseData> data, TaskType task_type)
{
    computeGlobal(data, DECAF_REDIST_SOURCE);

    splitData(data, DECAF_REDIST_SOURCE);

    redistribute(data, DECAF_REDIST_SOURCE);


}*/

void
decaf::
RedistZCurveMPI::computeGlobal(std::shared_ptr<BaseData> data, RedistRole role)
{
    if(role == DECAF_REDIST_SOURCE)
    {
        // If we don't have the global bounding box, we compute it once
        if(!bBBox_)
        {
            if(!data->hasZCurveKey())
            {
                std::cout<<"ERROR : Trying to redistribute the data with respect to a ZCurve "
                        <<"but no ZCurveKey (position) is available in the data. Abording."<<std::endl;
                MPI_Abort(MPI_COMM_WORLD, 0);
            }

            int nbParticules; // The size of the array is 3*nbParticules
            const float* key = data->getZCurveKey(&nbParticules);
            float localBBox[6];

            // Computing the local bounding box
            if(nbParticules > 0)
            {
                for(unsigned int i = 0; i < 6; i++)
                    localBBox[i] = key[i % 3];

                for(unsigned int i = 1; i < nbParticules; i++)
                {
                    for(unsigned int j = 0; j < 3; j++)
                    {
                        if(key[i*3+j] < localBBox[j])
                            localBBox[j] = key[i*3+j];
                        if(key[i*3+j] > localBBox[3+j])
                            localBBox[3+j] = key[i*3+j];
                    }
                }

                //Aggregating the bounding boxes to get the global one
                MPI_Allreduce(&(localBBox[0]), &(bBox_[0]), 3, MPI_FLOAT, MPI_MIN, commSources_);
                MPI_Allreduce(&(localBBox[0])+3, &(bBox_[0])+3, 3, MPI_FLOAT, MPI_MAX, commSources_);

                std::cout<<"The global bounding box is : ["<<localBBox[0]<<","<<localBBox[1]<<","<<localBBox[2]<<"]"
                        <<"["<<localBBox[3]<<","<<localBBox[4]<<","<<localBBox[5]<<"]"<<std::endl;

                bBBox_ = true;
            }
            else
            {
                std::cout<<"ERROR : no particules present on one process. Case not handled. Abording"<<std::endl;
                MPI_Abort(MPI_COMM_WORLD, 0);
            }
        }
    }
}

void
decaf::
RedistZCurveMPI::splitData(shared_ptr<BaseData> data, RedistRole role)
{
    std::cout<<"Spliting with the rank "<<rank_<<std::endl;
    if(role == DECAF_REDIST_SOURCE){

        //Compute the split vector and the destination ranks
        std::vector<std::vector<int> > split_ranges = std::vector<std::vector<int> >( nbDests_);
        int nbParticules;
        const float* pos = data->getZCurveKey(&nbParticules);

        // Create the array which represents where the current source will emit toward
        // the destinations rank. 0 is no send to that rank, 1 is send
        if( summerizeDest_) delete  summerizeDest_;
         summerizeDest_ = new int[ nbDests_];
        bzero( summerizeDest_,  nbDests_ * sizeof(int)); // First we don't send anything


        for(int i = 0; i < nbParticules; i++)
        {
            //Computing the cell of the particule
            int x = (pos[3*i] - bBox_[0]) / slicesDelta_[0];
            int y = (pos[3*i+1] - bBox_[1]) / slicesDelta_[1];
            int z = (pos[3*i+1] - bBox_[2]) / slicesDelta_[2];

            //Safety in case of wrong rounding
            if(x < 0) x = 0;
            if(y < 0) y = 0;
            if(z < 0) z = 0;
            if(x >= slices_[0]) x = slices_[0] - 1;
            if(y >= slices_[1]) y = slices_[1] - 1;
            if(z >= slices_[2]) z = slices_[2] - 1;

            unsigned int morton = Morton_3D_Encode_10bit(x,y,z);

            //Computing the destination rank
            int destination;
            if(morton < rankOffset_ * (indexes_per_dest_+1))
            {
                destination = morton / (indexes_per_dest_+1);
            }
            else
            {
                morton -= rankOffset_ * (indexes_per_dest_+1);
                destination = rankOffset_ + morton / indexes_per_dest_;
            }

            split_ranges[destination].push_back(i);

            //We won't send a message if we send to self
            if(destination + local_dest_rank_ != rank_)
                summerizeDest_[destination] = 1;

            destList_.push_back(destination + local_dest_rank_);

        }

        std::cout<<"Data will be split in "<<split_ranges.size()<<" chunks"<<std::endl;
        std::cout<<"Size of Destinaton list :  "<<destList_.size()<<" chunks"<<std::endl;

        splitChunks_ =  data->split( split_ranges );
        std::cout<<splitChunks_.size()<<" chunks have been produced"<<std::endl;

        std::cout<<"Serializing the chunks..."<<std::endl;
        for(unsigned int i = 0; i < splitChunks_.size(); i++)
        {
            // TODO : Check the rank for the destination.
            // Not necessary to serialize if overlapping
            if(!splitChunks_.at(i)->serialize())
                std::cout<<"ERROR : unable to serialize one object"<<std::endl;
        }

        // Everything is done, now we can clean the data.
        // Data might be rewriten if producers and consummers are overlapping
        data->purgeData();
    }
    else
        std::cout<<"Destination, nothing to do in the split"<<std::endl;
    std::cout<<"End of spliting"<<std::endl;

}

void
decaf::
RedistZCurveMPI::redistribute(std::shared_ptr<BaseData> data, RedistRole role)
{
    std::cout<<"Redistributing with the rank "<<rank_<<std::endl;
    // Sum the number of emission for each destination
    int *sum;
    if(role == DECAF_REDIST_SOURCE)
    {
        if(rank_ == local_source_rank_) sum = new int[ nbDests_];
        MPI_Reduce( summerizeDest_, sum,  nbDests_, MPI_INT, MPI_SUM,
                   local_source_rank_, commSources_);
    }
    std::cout<<"Reduce the sum done."<<std::endl;

    // Sending to the rank 0 of the destinations
    if(role == DECAF_REDIST_SOURCE && rank_ == local_source_rank_)
    {
        //for(int i = 0; i < nbDests_;i++)
        //    std::cout<<"Destination "<<i<<" : "<<sum[i]<<std::endl;
        MPI_Request req;
        reqs.push_back(req);
        MPI_Isend(sum,  nbDests_, MPI_INT,  local_dest_rank_, MPI_METADATA_TAG, communicator_,&reqs.back());

        std::cout<<"Sending the sum to the first destination ("<<rankDest_ -  rankSource_<<")"<<std::endl;
    }


    // Getting the accumulated buffer on the destination side
    int *destBuffer;
    if(role == DECAF_REDIST_DEST && rank_ ==  local_dest_rank_) //Root of destination
    {
        std::cout<<"Waiting for incomming array"<<std::endl;
        MPI_Status status;
        MPI_Probe(MPI_ANY_SOURCE, MPI_METADATA_TAG, communicator_, &status);
        if (status.MPI_TAG == MPI_METADATA_TAG)  // normal, non-null get
        {

            destBuffer = new int[ nbDests_];
            MPI_Recv(destBuffer,  nbDests_, MPI_INT, local_source_rank_, MPI_METADATA_TAG, communicator_, MPI_STATUS_IGNORE);
            std::cout<<"Receiving the sum  ("<<rankDest_ -  rankSource_<<")"<<std::endl;
        }

        //for(int i = 0; i < nbDests_;i++)
        //    std::cout<<"Reception : Destination "<<i<<" : "<<destBuffer[i]<<std::endl;
    }

    // Scattering the sum accross all the destinations
    int nbRecep;
    if(role == DECAF_REDIST_DEST)
    {
        std::cout<<"Scaterring with the rank "<<rank_<<std::endl;
        MPI_Scatter(destBuffer,  1, MPI_INT, &nbRecep, 1, MPI_INT, 0, commDests_);

        if(rank_ ==  local_dest_rank_) delete destBuffer;
        std::cout<<"Scattering done."<<std::endl;
    }

    // At this point, each source knows where they have to send data
    // and each destination knows how many message it should receive

    //Processing the data exchange
    if(role == DECAF_REDIST_SOURCE)
    {
        std::cout<<"Sending the data from the sources"<<std::endl;
        for(unsigned int i = 0; i <  destList_.size(); i++)
        {
            //Sending to self, we simply copy the string from the out to in
            if(destList_.at(i) == rank_)
            {
                transit = splitChunks_.at(i);
            }
            else
            {
                MPI_Request req;
                reqs.push_back(req);
                std::cout<<"Sending message of size : "<<splitChunks_.at(i)->getOutSerialBufferSize()
                        <<" to destination "<<destList_.at(i)<<std::endl;
                MPI_Isend( splitChunks_.at(i)->getOutSerialBuffer(),
                          splitChunks_.at(i)->getOutSerialBufferSize(),
                          MPI_BYTE, destList_.at(i), MPI_DATA_TAG, communicator_, &reqs.back());
            }
        }
        std::cout<<"End of sending messages"<<std::endl;
        // Cleaning the data here because synchronous send.
        // TODO :  move to flush when switching to asynchronous send
        splitChunks_.clear();
        destList_.clear();
    }

    if(role == DECAF_REDIST_DEST)
    {
        std::cout<<"Receiving the data."<<std::endl;
        for(int i = 0; i < nbRecep; i++)
        {
            MPI_Status status;
            MPI_Probe(MPI_ANY_SOURCE, MPI_DATA_TAG, communicator_, &status);
            if (status.MPI_TAG == MPI_DATA_TAG)  // normal, non-null get
            {
                int nitems; // number of items (of type dtype) in the message
                MPI_Get_count(&status, MPI_BYTE, &nitems);
                std::cout<<"Reception of a message with size "<<nitems<<std::endl;

                //Allocating the space necessary
                data->allocate_serial_buffer(nitems);
                //std::vector<char> buffer(nitems);
                std::cout<<"Allocation done"<<std::endl;
                MPI_Recv(data->getInSerialBuffer(), nitems, MPI_BYTE, status.MPI_SOURCE,
                         status.MPI_TAG, communicator_, &status);
                //MPI_Recv(&buffer[0], nitems, MPI_BYTE, status.MPI_SOURCE,
                //                         status.MPI_TAG, communicator_, &status);
                std::cout<<"Message received"<<std::endl;
                data->merge();
                //data->merge(&buffer[0], nitems);
                //A modifier
                //shared_ptr<BaseData> newData = make_shared<BaseData>(data->get());
                //receivedChunks_.push_back(newData);
            }

        }

        // Checking if we have something in transit
        if(transit)
        {
            std::cout<<"Getting the transit data"<<std::endl;
            data->merge(transit->getOutSerialBuffer(), transit->getOutSerialBufferSize());

            //We don't need it anymore. Cleaning for the next iteration
            transit.reset();
        }
        std::cout<<"End of reception."<<std::endl;
    }
    std::cout<<"End of redistributing with the rank "<<rank_<<std::endl;
}
// Merge the chunks from the vector receivedChunks into one single Data.
std::shared_ptr<decaf::BaseData>
decaf::
RedistZCurveMPI::merge(RedistRole role)
{
    return std::shared_ptr<BaseData>();
}




#endif