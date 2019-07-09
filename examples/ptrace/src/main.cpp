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

// #include "misc.h"



// consumer
void con(Decaf* decaf)
{   


		// mpaso mpas_c;
	vector< pConstructData > in_data;

	
	diy::mpi::communicator world(decaf->con_comm_handle());

	while (decaf->get(in_data))
	{
		int n_velocityX = 0;
		std::vector<int> indexToCellID;
		std::vector<double> data_bar;
		std::vector<int> data_bar_int;
		int fir_cell_idx = 999;
		// get the values and add them
		for (size_t i = 0; i < in_data.size(); i++)
		{


			SimpleFieldi d_metadata = in_data[i]->getFieldData<SimpleFieldi>("field_id");
			SimpleFieldi frame_num = in_data[i]->getFieldData<SimpleFieldi>("frame_num");

			// printf("fieldid %d \n", d_metadata.getData());


			switch(d_metadata.getData()){

				case 0: fprintf(stderr, "Recv xCell %d,\n", d_metadata.getData());
						break;

				case 1: fprintf(stderr, "Recv yCell %d,\n", d_metadata.getData());
						break;

				case 2:	fprintf(stderr, "Recv zCell %d,\n", d_metadata.getData());
						break;

				case 3:	fprintf(stderr, "Recv xVertex %d,\n", d_metadata.getData());
						break;
						
				case 4:	fprintf(stderr, "Recv yVertex %d,\n", d_metadata.getData());
						break;

				case 5:	fprintf(stderr, "Recv zVertex %d,\n", d_metadata.getData());
						break;

				case 6:	fprintf(stderr, "Recv indexToVertexID %d,\n", d_metadata.getData());
						break;

				case 7:	fprintf(stderr, "Recv indexToCellID %d,\n", d_metadata.getData());
						break;

				case 8:	fprintf(stderr, "Recv verticesOnEdge %d,\n", d_metadata.getData());
						break;

				case 9:	fprintf(stderr, "Recv cellsOnVertex %d,\n", d_metadata.getData());
						break;

				case 10:	fprintf(stderr, "Recv verticesOnCell %d,\n", d_metadata.getData());
						break;

				case 11:	fprintf(stderr, "Recv velocityXv %d, fno %d\n", d_metadata.getData(), frame_num.getData());
						break;

				case 12:	fprintf(stderr, "Recv velocityYv %d, fno %d\n", d_metadata.getData(), frame_num.getData());
						break;

				case 13:	fprintf(stderr, "Recv velocityZv %d,\n", d_metadata.getData());
						break;

				case 14:	fprintf(stderr, "Recv nEdgesOnCell %d,\n", d_metadata.getData());
						break;

				case 15:	fprintf(stderr, "Recv maxLevelCell %d,\n", d_metadata.getData());
						break;

				case 16:	fprintf(stderr, "Recv boundaryVertex %d,\n", d_metadata.getData());
						break;

				case 17:	fprintf(stderr, "Recv vertVelocityTop %d,\n", d_metadata.getData());
						break;

				case 18:	fprintf(stderr, "Recv cellsOnCell %d,\n", d_metadata.getData());
						break;

				case 19:	fprintf(stderr, "Recv zTop %d,\n", d_metadata.getData());
						break;

				case 20:	fprintf(stderr, "Recv zMid %d,\n", d_metadata.getData());
						break;

				default:fprintf(stderr, "Exiting! %d,\n", d_metadata.getData()); 
						exit(0);
						break;


			}

			// if(d_metadata.getData()==19)
			// {	
			// 	VectorFliedd d_data_bar = in_data[i]->getFieldData<VectorFliedd>("data_d");

			// 	data_bar = d_data_bar.getVector();
			// 	printf("Incomingg container %f metadata2 %d\n", data_bar[data_bar.size()-1], d_metadata.getData());

				
			// }else if (d_metadata.getData()==7){
			// 	VectorFieldi d_data_bar = in_data[i]->getFieldData<VectorFieldi>("data_i");

			// 	data_bar_int = d_data_bar.getVector();
			// 	printf("Incomingg container %d metadata2 %d\n", data_bar_int[data_bar.size()-1], d_metadata.getData());

			// }else if (d_metadata.getData()==0){

			// 	VectorFliedd d_data_bar = in_data[i]->getFieldData<VectorFliedd>("data_d");
			// 	data_bar = d_data_bar.getVector();
			// 	printf("Incomingg container %f metadata2 %d\n", data_bar[data_bar.size()-1], d_metadata.getData());
			// }
			

			
		}
	}

}

int main(){

		// define the workflow
	Workflow workflow;
	//make_wflow(workflow);
	Workflow::make_wflow_from_json(workflow, "mpas_decaf_flowvis.json");

	MPI_Init(NULL, NULL);


	// create decaf
	Decaf* decaf = new Decaf(MPI_COMM_WORLD, workflow);

	con(decaf);


	delete decaf;

	MPI_Finalize();

}