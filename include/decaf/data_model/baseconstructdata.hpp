#ifndef BASE_CONSTRUCT_DATA
#define BASE_CONSTRUCT_DATA

#include <decaf/decaf.hpp>
#include <decaf/data_model/basedata.h>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include "serialize_tuple.h"
//#include <boost/serialization/shared_ptr.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/unordered_map.hpp>



#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <map>
#include <vector>
#include <sstream>

enum ConstructTypeFlag {
    DECAF_NOFLAG = 0x0,     // No specific information on the data field
    DECAF_NBITEM = 0x1,     // Field represents the number of item in the collection
    DECAF_ZCURVEKEY = 0x2,  // Field that can be used as a key for the ZCurve (position)
    DECAF_ZCURVEINDEX = 0x4 // Field that can be used as the index for the ZCurve (hilbert code)
};

enum ConstructTypeScope {
    DECAF_SHARED = 0x0,     // This value is the same for all the items in the collection
    DECAF_PRIVATE = 0x1,    // Different values for each items in the collection
};

enum ConstructTypeSplitPolicy {
    DECAF_SPLIT_DEFAULT = 0x0,      // Call the split fonction of the data object
    DECAF_SPLIT_KEEP_VALUE = 0x1,   // Keep the same values for each split
    DECAF_SPLIT_MINUS_NBITEM = 0x2, // Withdraw the number of items to the current value
};

enum ConstructTypeMergePolicy {
    DECAF_MERGE_DEFAULT = 0x0,        // Call the split fonction of the data object
    DECAF_MERGE_FIRST_VALUE = 0x1,    // Keep the same values for each split
    DECAF_MERGE_ADD_VALUE = 0x2,      // Add the values
    DECAF_MERGE_APPEND_VALUES = 0x4,  // Append the values into the current object
};

//---/ Wrapper for std::shared_ptr<> /------------------------------------------

namespace boost { namespace serialization {

template<class Archive, class Type>
void save(Archive & archive, const std::shared_ptr<Type> & value, const unsigned int /*version*/)
{
    Type *data = value.get();
    archive << data;
}

template<class Archive, class Type>
void load(Archive & archive, std::shared_ptr<Type> & value, const unsigned int /*version*/)
{
    Type *data;
    archive >> data;

    typedef std::weak_ptr<Type> WeakPtr;
    static boost::unordered_map<void*, WeakPtr> hash;

    if (hash[data].expired())
    {
         value = std::shared_ptr<Type>(data);
         hash[data] = value;
    }
    else value = hash[data].lock();
}

template<class Archive, class Type>
inline void serialize(Archive & archive, std::shared_ptr<Type> & value, const unsigned int version)
{
    split_free(archive, value, version);
}

}}

class BaseConstructData {
public:
    BaseConstructData(){}
    virtual ~BaseConstructData(){}

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    }

    virtual int getNbItems() = 0;

    virtual std::vector<std::shared_ptr<BaseConstructData> > split(
            const std::vector<int>& range,
            ConstructTypeSplitPolicy policy = DECAF_SPLIT_DEFAULT ) = 0;

    virtual std::vector<std::shared_ptr<BaseConstructData> > split(
            const std::vector< std::vector<int> >& range,
            ConstructTypeSplitPolicy policy = DECAF_SPLIT_DEFAULT ) = 0;

    virtual bool merge(std::shared_ptr<BaseConstructData> other,
                       ConstructTypeMergePolicy policy = DECAF_MERGE_DEFAULT) = 0;

    virtual bool canMerge(std::shared_ptr<BaseConstructData> other) = 0;

};
BOOST_CLASS_EXPORT_GUID(BaseConstructData,"BaseConstructData")

#endif