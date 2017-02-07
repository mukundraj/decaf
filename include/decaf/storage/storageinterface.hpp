//---------------------------------------------------------------------------
// Define the interface of a storage object to store frames
//
// Matthieu Dreher
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// mdreher@anl.gov
//
//--------------------------------------------------------------------------

#ifndef DECAF_STORAGE_INTERFACE
#define DECAF_STORAGE_INTERFACE

#include <decaf/data_model/pconstructtype.h>

namespace decaf
{
    class Storage
    {
    public:
        Storage(){}

        virtual ~Storage(){}

        virtual bool isFull() = 0;
        virtual unsigned int getBufferSize() = 0;
        virtual bool insert(unsigned int id, pConstructData data) = 0;
        virtual void erase(unsigned int id) = 0;
        virtual bool hasData(unsigned int id) = 0;
        virtual pConstructData getData(unsigned int id) = 0;
    };
} // namespace decaf

#endif
