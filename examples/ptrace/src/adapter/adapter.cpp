#include "adapter.h"

#include <decaf/decaf.hpp>
#include <bredala/data_model/pconstructtype.h>
#include <bredala/data_model/simplefield.hpp>
#include <bredala/data_model/arrayfield.hpp>
#include <bredala/data_model/vectorfield.hpp>
#ifdef TRANSPORT_CCI
#include <cci.h>
#endif
#include <bredala/data_model/boost_macros.h>

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <map>
#include <cstdlib>


using namespace decaf;
using namespace std;


Decaf* decaf2;


void f_MPI_Comm_c2f(MPI_Fint &comm) {

	comm  = decaf2->prod_comm_handle();

}

void decaf_put_int_array(int id, int frame_no, int size, int i1d[]){

	pConstructData container ;

	SimpleFieldi field_id(id);
	container->appendData("field_id", field_id, DECAF_NOFLAG, DECAF_PRIVATE, DECAF_SPLIT_DEFAULT, DECAF_MERGE_DEFAULT); 

	VectorFieldi data (&i1d[0], size, size);
	container->appendData("data_i", data, DECAF_NOFLAG, DECAF_PRIVATE, DECAF_SPLIT_DEFAULT, DECAF_MERGE_DEFAULT); 

	double garbage;
	VectorFliedd data_d (&garbage, 1, 1);
	container->appendData("data_d", data_d, DECAF_NOFLAG, DECAF_PRIVATE, DECAF_SPLIT_DEFAULT, DECAF_MERGE_DEFAULT);

	decaf2->put(container);

}


void decaf_put_double_array(int id, int frame_no, int size, double d1d[]){

	pConstructData container ;



	SimpleFieldi field_id(id);
	container->appendData("field_id", field_id, DECAF_NOFLAG, DECAF_PRIVATE, DECAF_SPLIT_DEFAULT, DECAF_MERGE_DEFAULT); 

	SimpleFieldi frame_num(frame_no);
	container->appendData("frame_num", frame_num, DECAF_NOFLAG, DECAF_PRIVATE, DECAF_SPLIT_DEFAULT, DECAF_MERGE_DEFAULT); 

	VectorFliedd data (&d1d[0], size, size);
	container->appendData("data_d", data, DECAF_NOFLAG, DECAF_PRIVATE, DECAF_SPLIT_DEFAULT, DECAF_MERGE_DEFAULT); 

	int garbage;
	VectorFieldi data_i (&garbage, 1, 1);
	container->appendData("data_i", data_i, DECAF_NOFLAG, DECAF_PRIVATE, DECAF_SPLIT_DEFAULT, DECAF_MERGE_DEFAULT);

	decaf2->put(container);

}

void init_decaf(){

	
	Workflow workflow;
	Workflow::make_wflow_from_json(workflow, "mpas_decaf_flowvis.json");
	MPI_Init(NULL, NULL);
	decaf2 = new Decaf(MPI_COMM_WORLD, workflow);
	//MPI_Fint fcommunicator = MPI_Comm_c2f(decaf2->prod_comm_handle());
	//communicator = &fcommunicator;
	// run decaf
	//run(workflow, MPI_CO);

	//printf("pid %d initialized \n", decaf2->world->rank());
}



void finish_decaf(){

	if (decaf2->my_node("prod") || decaf2->my_node("con")){
		int pid=decaf2->world->rank();
		printf("pid %d terminating now \n", decaf2->world->rank());
		decaf2->terminate();
	}

	// int flag = false;
	// MPI_Finalized(&flag);
	// if (!flag)
	// MPI_Finalize();
}


// link callback function
extern "C"
{
	// dataflow just forwards everything that comes its way in this example
	void dflow(void* args,                          // arguments to the callback
			Dataflow* dataflow,                  // dataflow
			pConstructData in_data)   // input data
	{
		int var = in_data->getFieldData<SimpleFieldi>("var").getData();
		fprintf(stderr, "Forwarding data %d in dflow\n", var);
		dataflow->put(in_data, DECAF_LINK);
	}
} // extern "C"

