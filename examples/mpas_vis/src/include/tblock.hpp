

//
// block structure
// the contents of a block are completely user-defined
//
struct TBlock
{
    std::vector<int>      values;
    float                 average;

};

//
// When a block has mutliple fields of different types -- specifically, when it
// cannot be serialized by simple binary copying -- diy needs a little more
// information about how to serialize the block. This is done by specializing
// diy::Serialization class for the block type and defining the appropriate
// save and load functions. Usually, this amounts to listing the member
// variables; diy knows how to serialize simple (POD) types and some STL types.
//
namespace diy
{
    template<>
        struct Serialization<TBlock>
	{
        static void save(BinaryBuffer& bb, const TBlock& b)
	{
            diy::save(bb, b.values);
            diy::save(bb, b.average);
        
	}

        static void load(BinaryBuffer& bb, TBlock& b)
	{
            diy::load(bb, b.values);
            diy::load(bb, b.average);
        
	}
    
	};

}

//
// Functions must be defined to create, destroy, save, and load a block.
// Create and destroy allocate and free the block, while
// save and load serialize and deserialize the block.
// These four functions are called when blocks are cycled in- and out-of-core.
// They could have been static member functions of the block above, or kept
// separate, as below
//
void* create_block()
{
    TBlock* b = new TBlock;
    return b;

}

void destroy_block(void* b)
{
    delete static_cast<TBlock*>(b);

}

void save_block(const void* b, diy::BinaryBuffer& bb)
{
    diy::save(bb, *static_cast<const TBlock*>(b));

}

void load_block(void* b, diy::BinaryBuffer& bb)
{
    diy::load(bb, *static_cast<TBlock*>(b));

}
