
#include <decaf/decaf.hpp>
#include <bredala/data_model/simplefield.hpp>
#include <bredala/data_model/vectorfield.hpp>
#include <bredala/data_model/boost_macros.h>

#include <diy/mpi.hpp>
#include <diy/master.hpp>
#include <diy/reduce-operations.hpp>
#include <diy/decomposition.hpp>
#include <diy/mpi/datatypes.hpp>
#include <diy/io/bov.hpp>
#include <diy/pick.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/merge.hpp>

#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <map>
#include <cstdlib>

using namespace decaf;
using namespace std;



#include<iostream>
#include<string>
#include "streamlines.h"
#include "pathlines.h"
#include "mpaso.h" 

#include "ptrace.hpp"
#include "block.hpp"
#include "tblock.hpp"

int frame_number = 0;

// add a block to the master
struct AddBlock
{
    AddBlock(diy::Master &master_) :
        master(master_)
    {}

    DBlock* operator()(int gid,
                      const diy::Link& link) const
    {
        DBlock *b       = new DBlock;
        diy::Link *l      = new diy::Link(link);
        diy::Master &m = const_cast<diy::Master&>(master);
        m.add(gid, b, l);
        return b;
    
    }

    diy::Master& master;

};


// add a block to the master and read input data
struct AddAndRead : public AddBlock
{
    AddAndRead(diy::Master& m,
               const char*  infile_,
               diy::mpi::communicator& world_,
               const float vec_scale_,
               const int hdr_bytes_) :
        AddBlock(m),
        infile(infile_),
        world(world_),
        vec_scale(vec_scale_),
        hdr_bytes(hdr_bytes_) {}

    DBlock* operator()(int gid,
                    const diy::Link& link) const
    {
        DBlock* b = AddBlock::operator()(gid, link);
	return b;
        /*MPI_Offset *start, *count;
        float *data_u=NULL, *data_v=NULL, *data_w=NULL;

        int ncfile, ndims, nvars, ngatts, unlimited;
        int ret;
        ret = ncmpi_open(world, infile, NC_NOWRITE, MPI_INFO_NULL,&ncfile);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        ret = ncmpi_inq(ncfile, &ndims, &nvars, &ngatts, &unlimited);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);


        // reversed order of shape and bounds needed because the sample data file
        // is linearized in row-major (C) order
        vector<int> shape(3);
        for (size_t i = 0; i < 3; i++)
            shape[2 - i] = domain.max[i] - domain.min[i] + 1;
        //        diy::io::BOV reader(in, shape, hdr_bytes);

        Bounds r_bounds;
        r_bounds.min[0] = bounds.min[2];
        r_bounds.max[0] = bounds.max[2];
        r_bounds.min[1] = bounds.min[1];
        r_bounds.max[1] = bounds.max[1];
        r_bounds.min[2] = bounds.min[0];
        r_bounds.max[2] = bounds.max[0];



        start = (MPI_Offset*) calloc(ndims, sizeof(MPI_Offset));
        count = (MPI_Offset*) calloc(ndims, sizeof(MPI_Offset));

	if (ndims==4){
            count[0] = 1;
            count[1] = r_bounds.max[0] - r_bounds.min[0]+1;
            count[2] = r_bounds.max[1] - r_bounds.min[1]+1;
            count[3] = r_bounds.max[2] - r_bounds.min[2]+1;

            start[0] =  0; start[1] = r_bounds.min[0]; start[2] = r_bounds.min[1]; start[3] = r_bounds.min[2];
        
	}else if(ndims==3){

            count[0] = r_bounds.max[0] - r_bounds.min[0]+1;
            count[1] = r_bounds.max[1] - r_bounds.min[1]+1;
            count[2] = r_bounds.max[2] - r_bounds.min[2]+1;

            start[0] = r_bounds.min[0]; start[1] = r_bounds.min[1]; start[2] = r_bounds.min[2];
        
	}

        //        std::cout<<"counts"<<count[0]<<" "<<count[1]<<" "<<count[2]<<"\n";
        //        std::cout<<"starts"<<start[0]<<" "<<start[1]<<" "<<start[2]<<"\n";

        size_t nvecs =
                (bounds.max[0] - bounds.min[0] + 1) *
                (bounds.max[1] - bounds.min[1] + 1) *
                (bounds.max[2] - bounds.min[2] + 1);
        vector<float> values(nvecs * 3); // temporary contiguous buffer of input vector values

        data_u = (float*) calloc(nvecs, sizeof(float));
        data_v = (float*) calloc(nvecs, sizeof(float));
        data_w = (float*) calloc(nvecs, sizeof(float));
        ret = ncmpi_get_vara_float_all(ncfile, 0, start, count, data_u);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
        ret = ncmpi_get_vara_float_all(ncfile, 1, start, count, data_v);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
        ret = ncmpi_get_vara_float_all(ncfile, 2, start, count, data_w);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);



        // copy from temp values into block
        b->vel[0] = new float[nvecs];
        b->vel[1] = new float[nvecs];
        b->vel[2] = new float[nvecs];
        b->nvecs = nvecs;
        for (size_t i = 0; i < nvecs; i++)
	{

            b->vel[0][i] = data_u[i] * vec_scale;
            b->vel[1][i] = data_v[i] * vec_scale;
            b->vel[2][i] = data_w[i] * vec_scale;

        
	}

        ret = ncmpi_close(ncfile);
        free(start);
        free(count);
        free(data_u);
        free(data_v);
        free(data_w);

*/
    
    }


    const char* infile;
    diy::mpi::communicator world;
    float vec_scale;
    int hdr_bytes;

};


// consumer
void con(Decaf* decaf)
{       
	mpaso mpas1;
	vector< pConstructData > in_data;


	diy::mpi::communicator world(decaf->con_comm_handle());
	// defaults
	int nblocks     = world.size();           // total number of global blocks
	int nthreads    = 1;                      // number of threads diy can use
	int mblocks     = -1;                     // number of blocks in memory (-1 = all)
	string prefix   = "./DIY.XXXXXX";         // storage of temp files
	int ndims       = 3;                      // domain dimensions
	float vec_scale = 1.0;                    // vector field scaling factor
	int hdr_bytes   = 0;                      // num bytes header before start of data in infile
	int max_rounds  = 0;                      // max number of rounds to trace (0 = no limit)

	diy::FileStorage             storage(prefix);	
	diy::Master                     master(world,
			nthreads,
			mblocks,
			&DBlock::create,
			&DBlock::destroy,
			&storage,
			&DBlock::save,
			&DBlock::load);
	diy::RoundRobinAssigner     assigner(world.size(), nblocks);

	string infile;
	AddAndRead addblock(master, infile.c_str(), world, vec_scale, hdr_bytes);	
	
	std::vector<int> gids;                     // global ids of local blocks
	assigner.local_gids(world.rank(), gids);   // get the gids of local blocks
	for (unsigned i = 0; i < gids.size(); ++i) // for the local blocks in this processor
	{       		
		int gid = gids[i];
		fprintf(stderr, " inside consumer, gid: %d \n", gid);

		//diy::Link    link;  // link is this block's neighborhood
	
		diy::Link* link = new diy::Link;
		/*
		diy::BlockID  neighbor;
		if (gid < nblocks - 1)               // all but the last block in the global domain
		{
		neighbor.gid  = gid + 1;                     // gid of the neighbor block
		neighbor.proc = assigner.rank(neighbor.gid); // process of the neighbor block
			link->add_neighbor(neighbor);                // add the neighbor block to the link

		}
		if (gid > 0)                         // all but the first block in the global domain
		{
			neighbor.gid  = gid - 1;
			neighbor.proc = assigner.rank(neighbor.gid);
			link->add_neighbor(neighbor);

		}
		
		Block* b = new Block;                // create a new block
		// assign some values for the block
		// in this example, simply based on the block global id
		for (unsigned j = 0; j < 3; ++j)
			b->values.push_back(gid * 3 + j);
		*/
		

		//DBlock* b = addblock(gid, *link);
		TBlock *b = new TBlock;
		master.add(gid, b, link);            // add the current local block to the master

	}



	//fprintf(stderr, "diy world size: %d", world.size() );
	std::string ip_file = "graph.topology";
	mpas1.generate_domain_decomposition_graph(ip_file, world.size());	

	int ctr = 0;
	while (decaf->get(in_data))
	{
		int n_velocityX = 0;
		std::vector<int> indexToCellID;
		std::vector<double> data_bar;
		int fir_cell_idx = 999;
		// get the values and add them
		for (size_t i = 0; i < in_data.size(); i++)
		{

			VectorFliedd d_data_bar = in_data[i]->getFieldData<VectorFliedd>("data_bar");
			if (d_data_bar){
				data_bar = d_data_bar.getVector();
				printf("got data_bar in con: %f %f %f %f %f %f %f %f %f %f\n", data_bar[13], data_bar[14], data_bar[15], data_bar[16], data_bar[17], data_bar[18], data_bar[1452848+0], data_bar[1452848+4], data_bar[1452848+5], data_bar[1452848+6]);
				fir_cell_idx = data_bar[0];
				//printf("fir cell val %f\n", data_bar[int(12+data_bar[4]+data_bar[5]+data_bar[6]+data_bar[7])]);
				fprintf(stderr, "databar size :%ld\n", data_bar.size());
				unsigned int init_t=0, fin_t = 40000000, h = 200000;
				std::vector<double> seeds_xyz(1000*3);
				std::string op_file = std::to_string(ctr)+".vtp";
				//streamlines slines(mpas1, data_bar, op_file , init_t, fin_t, h, seeds_xyz);
			}else{ 
				printf("null ptr data_bar\n");
			}
		}
		// printf("consumer sum = %d\n", sum);
		ctr++;
	}


	// terminate the task (mandatory) by sending a quit message to the rest of the workflow
	fprintf(stderr, "consumer terminating\n");
	decaf->terminate();

}

int main(){


	std::cout<<"starting mpas consumer\n";

	// define the workflow
	Workflow workflow;
	//make_wflow(workflow);
	Workflow::make_wflow_from_json(workflow, "mpas_decaf_flowvis.json");

	MPI_Init(NULL, NULL);
	//diy::mpi::environment env(NULL, NULL);

	// create decaf
	Decaf* decaf = new Decaf(MPI_COMM_WORLD, workflow);

        
	// start the task
	con(decaf);

	// cleanup
	delete decaf;

	fprintf(stderr, "post delete decaf\n");
	MPI_Finalize();

	/*
	   std::cout.precision(17);

	   int nseeds = 1;
	//    std::string ip_file = "/Users/mukundraj/Desktop/work/datasets/output.0014-05-01_00.00.00.nc";
	std::string ip_file = "/homes/mraj/work/projects/MPAS-Model-6.0/testdir/MPAS-O_V6.0_QU240/output.nc";
	std::string op_file = "";
	//    unsigned int init_t=0, fin_t = 691199, h = 10000, interval = 172800;
	unsigned int init_t=0, fin_t = 40000000, h = 200000; // iters = fin_t/h
	std::vector<double> seeds_xyz(nseeds*3);

	streamlines slines(ip_file, , init_t, fin_t, h, seeds_xyz);
	//    pathlines plines(ip_file, op_file, init_t, fin_t, interval, h, seeds_xyz);

	 */



	std::cout<<"finished2\n";

	return 0;

}
