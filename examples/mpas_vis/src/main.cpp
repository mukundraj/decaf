
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

int frame_number = 0;



// consumer
void con(Decaf* decaf)
{       
	mpaso mpas1;
	vector< pConstructData > in_data;
	

	diy::mpi::communicator world(decaf->con_comm_handle());

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
				//printf("got data_bar in con: %f %f %f %f %f %f %f %f %f %f\n", data_bar[13], data_bar[14], data_bar[15], data_bar[16], data_bar[17], data_bar[18], data_bar[1452848+0], data_bar[1452848+4], data_bar[1452848+5], data_bar[1452848+6]);
				fir_cell_idx = data_bar[0];
				//printf("fir cell val %f\n", data_bar[int(12+data_bar[4]+data_bar[5]+data_bar[6]+data_bar[7])]);
				fprintf(stderr, "databar size :%ld\n", data_bar.size());
				unsigned int init_t=0, fin_t = 40000000, h = 200000;
				std::vector<double> seeds_xyz(1000*3);
				std::string op_file = std::to_string(ctr)+".vtp";
				//streamlines slines(mpas1, data_bar, op_file , init_t, fin_t, h, seeds_xyz);
			}else{ 
				fprintf(stderr, "null ptr data_bar\n");
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

	// create decaf
	Decaf* decaf = new Decaf(MPI_COMM_WORLD, workflow);

	// start the task
	con(decaf);

	// cleanup
	delete decaf;
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



	std::cout<<"finished\n";

	return 0;

}
