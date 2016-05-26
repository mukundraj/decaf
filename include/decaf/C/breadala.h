#ifndef BREDALA_H
#define BREDALA_H

/* For compatibility with Clang. */
#ifndef __has_extension
#	define __has_extension(EXT) 0
#endif

#ifdef __cplusplus
#include <cstddef>
using std::size_t;
extern "C" {
#else /* __cplusplus */
#include <stdbool.h>
#include <stddef.h>
#endif /* __cplusplus */


    /* For each  type a pointer  type to an  unspecified (in fact  not existing)
     * struct is created. In fact, they are pointers to C++ objects.
     *
     * These intermediary structs are used to enforce type safety. */
#   define BCA_DECLARE_TYPE(NAME) \
    typedef struct NAME ## _inexistent *NAME

    /**
     * A constructdata is a container to store all the fields of a data model
     */
    BCA_DECLARE_TYPE(bca_constructdata);

    /**
     * A simple field is a straight object which can not be subdivided further by
     * the system. (ex : int, float, etc)
     * The field should be pushed in a bca_constructdata
     */
    //BCA_DECLARE_TYPE(bca_simplefield);

    /**
     * Field type to store an array. The array may be subdivided by the system but
     * not the individual elements.
     * The field should be pushed in a bca_constructdata
     */
    //BCA_DECLARE_TYPE(bca_arrayfield);

    BCA_DECLARE_TYPE(bca_field);


    /**
     * Flag to define the type of simplefield or arrayfield.
     */
    typedef enum {
         bca_INT,
         bca_UNSIGNED,
         bca_FLOAT,
         bca_DOUBLE,
         bca_CHAR
    } bca_type;

    // Has to be exactly as in maptools.h
    typedef enum {
        bca_NOFLAG = 0x0,     // No specific information on the data field
        bca_NBITEM = 0x1,     // Field represents the number of item in the collection
        bca_ZCURVEKEY = 0x2,  // Field that can be used as a key for the ZCurve (position)
        bca_ZCURVEINDEX = 0x4 // Field that can be used as the index for the ZCurve (hilbert code)
    } bca_ConstructTypeFlag;

    typedef enum  {
        bca_SHARED = 0x0,     // This value is the same for all the items in the collection
        bca_PRIVATE = 0x1,    // Different values for each items in the collection
        bca_SYSTEM = 0x2,     // This field is not a user data
    } bca_ConstructTypeScope;

    typedef enum {
        bca_SPLIT_DEFAULT = 0x0,      // Call the split fonction of the data object
        bca_SPLIT_KEEP_VALUE = 0x1,   // Keep the same values for each split
        bca_SPLIT_MINUS_NBITEM = 0x2, // Withdraw the number of items to the current value
        bca_SPLIT_SEGMENTED = 0x4,
    } bca_ConstructTypeSplitPolicy;

    typedef enum  {
        bca_MERGE_DEFAULT = 0x0,        // Call the split fonction of the data object
        bca_MERGE_FIRST_VALUE = 0x1,    // Keep the same values for each split
        bca_MERGE_ADD_VALUE = 0x2,      // Add the values
        bca_MERGE_APPEND_VALUES = 0x4,  // Append the values into the current object
        bca_MERGE_BBOX_POS = 0x8,       // Compute the bounding box from the field pos
    } bca_ConstructTypeMergePolicy;

    bca_constructdata
    bca_new_constructdata();

    bool
    bca_append_field(   bca_constructdata container,
                        const char* name,
                        bca_field field,
                        bca_ConstructTypeFlag flag,
                        bca_ConstructTypeScope scope,
                        bca_ConstructTypeSplitPolicy split,
                        bca_ConstructTypeMergePolicy merge);

    bca_field
    bca_create_simplefield(void* data,
                           bca_type type);

    bca_field
    bca_create_arrayfield(void* data,
                          bca_type type,
                          int size,
                          int element_per_items,
                          int capacity,
                          bool owner);

    void bca_free_constructdata(bca_constructdata data);

    void bca_free_field(bca_field data);


    bca_field
    bca_get_simplefield(bca_constructdata container, const char* name, bca_type type);

    bca_field
    bca_get_arrayfield(bca_constructdata container, const char* name, bca_type type);

    bool
    bca_get_data(bca_field field, bca_type type, void* data);

    void * bca_get_array(bca_field field, bca_type type, size_t* size);

    int
    bca_get_nbitems(bca_field field);




#undef BCA_DEPRECATED

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* BREDALA_H */

