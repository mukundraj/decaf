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

int *gar, num_arrived_fields = 0;


Decaf* decaf2;
pConstructData container;


double *velocityX, *velocityY, *velocityZ, *zTop;
int n_velocityX, n_velocityY, n_velocityZ, n_zTop;
int *indexToCellID, n_indexToCellID, *indexToVertexID, n_indexToVertexID, *cellsOnVertex, n_cellsOnVertex;
std::vector<double> db_indexToCellID, db_indexToVertexID, db_cellsOnVertex;
double *xCell, *yCell, *zCell, *xVertex, *yVertex, *zVertex;
int n_xCell, n_yCell, n_zCell, n_xVertex, n_yVertex, n_zVertex;


void func1(int n, double x[]){

	int i=0;

}

void decaf_put(int id, int frame_no, int nVertices, int nCells, int nVertLevels, int n_procs){
	//printf("decaf_put called by fieldid: %d\n", id);
	pConstructData container ;

	std::vector<double> data_bar(20);
	int n_data_bar = 20;


	//VectorFliedd data3 (velocityX, n_velocityX, 1);
	//container->appendData("velocityX", data3, DECAF_NOFLAG, DECAF_PRIVATE, DECAF_SPLIT_DEFAULT, DECAF_MERGE_DEFAULT); 

	//VectorFliedd data2 (xCell, n_xCell, 1);
	//container2->appendData("xCell", data2, DECAF_NOFLAG, DECAF_PRIVATE, DECAF_SPLIT_DEFAULT, DECAF_MERGE_DEFAULT);

	
	//VectorFieldi data1 (indexToCellID, n_indexToCellID, 1);
	// container->appendData("indexToCellID", data1, DECAF_NOFLAG, DECAF_PRIVATE, DECAF_SPLIT_DEFAULT, DECAF_MERGE_DEFAULT);

	// velocityX 4
	data_bar.insert(data_bar.end(), &velocityX[0], &velocityX[n_velocityX]);
	data_bar[4] = double(n_velocityX);
	n_data_bar += n_velocityX;

	// velocityY 5
	data_bar.insert(data_bar.end(), &velocityY[0], &velocityY[n_velocityY]);
	data_bar[5] = double(n_velocityY);
	n_data_bar += n_velocityY;

	// velocityZ 6
	data_bar.insert(data_bar.end(), &velocityZ[0], &velocityZ[n_velocityZ]);
	data_bar[6] = double(n_velocityZ);
	n_data_bar += n_velocityZ;

	// zTop 11
	data_bar.insert(data_bar.end(), &zTop[0], &zTop[n_zTop]);
	data_bar[11] = double(n_zTop);
	n_data_bar += n_zTop;

	// indexToCellID 0
	data_bar.insert(data_bar.end(), db_indexToCellID.begin(), db_indexToCellID.end());
	data_bar[0] = double(n_indexToCellID);
	n_data_bar += n_indexToCellID;
      

 
	//printf("data_bar first five vals + size  %d %d %d %d %d %d %d\n", int(data_bar[0]), int(data_bar[1]), int(data_bar[2]), int(data_bar[3]), int(data_bar[4]), int(data_bar.size()), n_data_bar);
 
        

	if (id==12){
		// xCell 1
		//printf("size before xCell %d. xCell[0, 1, 2] %f %f %f\n", int(data_bar.size()), xCell[0], xCell[1], xCell[2]);
		data_bar.insert(data_bar.end(), &xCell[0], &xCell[n_xCell]);
		data_bar[1] = double(n_xCell);
		n_data_bar += n_xCell;

		// yCell 2
		data_bar.insert(data_bar.end(), &yCell[0], &yCell[n_yCell]);
		data_bar[2] = double(n_yCell);
		n_data_bar += n_yCell;

		// zCell 3
		data_bar.insert(data_bar.end(), &zCell[0], &zCell[n_zCell]);
		data_bar[3] = double(n_zCell);
		n_data_bar += n_zCell;

		// xVertex 8
		data_bar.insert(data_bar.end(), &xVertex[0], &xVertex[n_xVertex]);
		data_bar[8] = double(n_xVertex);
		n_data_bar += n_xVertex;

		// yVertex 9
		data_bar.insert(data_bar.end(), &yVertex[0], &yVertex[n_yVertex]);
		data_bar[9] = double(n_yVertex);
		n_data_bar += n_xVertex;

		// zVertex 10
		data_bar.insert(data_bar.end(), &zVertex[0], &zVertex[n_zVertex]);
		data_bar[10] = double(n_zVertex);
		n_data_bar += n_zVertex;

	        // indexToVertexID 7
		data_bar.insert(data_bar.end(), db_indexToVertexID.begin(), db_indexToVertexID.end());
		data_bar[7] = double(n_indexToVertexID);
		n_data_bar += n_indexToVertexID;


		// cellsOnVertex 12
		data_bar.insert(data_bar.end(), db_cellsOnVertex.begin(), db_cellsOnVertex.end());
		data_bar[12] = double(n_cellsOnVertex);
		n_data_bar += n_cellsOnVertex;

	}
	

	data_bar[13] = double (frame_no);
	data_bar[14] = double (nVertices);
	data_bar[15] = double (nCells);
	data_bar[16] = double (nVertLevels);
	data_bar[17] = double (n_procs);
	data_bar[18] = double (data_bar.size());
	//fprintf(stderr, "sending bar %f \n", data_bar[0]);
	//printf("first xCell %f %f %f\n", data_bar[1445620], data_bar[1445621], data_bar[1445622]);
	VectorFliedd data (&data_bar[0], n_data_bar, n_data_bar); 
	container->appendData("data_bar", data, DECAF_NOFLAG, DECAF_PRIVATE, DECAF_SPLIT_DEFAULT, DECAF_MERGE_DEFAULT); 
	decaf2->put(container);
	
	// clear vectors
	//db_indexToVertexID.clear();
	db_indexToCellID.clear();
	data_bar.clear();
}

void stage_put_intarray(int i1d[], int n, int id, int frame_num){

	//            ArrayFieldi data(i1d, n, 1, true);
	//printf("my len: %d\n ", n);
	//pConstructData container;
	switch(id){

		case 0: indexToCellID = &i1d[0];
			//printf("indexToCellID [0]: %d\n", indexToCellID[0]);
			n_indexToCellID = n;
			db_indexToCellID.insert(db_indexToCellID.end(), &indexToCellID[0], &indexToCellID[n]);
			// fprintf(stderr, "indexToCell adapter %d \n", indexToCellID[0]);
			//container->appendData("indexToCellID", data, DECAF_NOFLAG, DECAF_PRIVATE, DECAF_SPLIT_DEFAULT, DECAF_MERGE_DEFAULT);
			break;
		case 7: indexToVertexID = &i1d[0];
			db_indexToVertexID.insert(db_indexToVertexID.end(), &indexToVertexID[0], &indexToVertexID[n]);
			n_indexToVertexID = n;
			break;          

		case 12: cellsOnVertex = &i1d[0];
			 n_cellsOnVertex = n;
			 db_cellsOnVertex.insert(db_cellsOnVertex.end(), &cellsOnVertex[0], &cellsOnVertex[n]);
			 //printf("cellsOnVertex ");
			//for (int i=0; i<10; i++){
			//	printf("%d ", cellsOnVertex[i]);
			//}

			break; 
		default:
			 break;     
	}

	//   if (num_arrived_fields == 2){
	//      decaf2->put(container);
	//	num_arrived_fields = 0;
	//   }
}
void stage_put_doublearray(double d1d[], int n, int id, int frame_num){
	// printf ("in stage_put_double\n");
	//printf("my len d: %d\n ", n);
	//pConstructData container;
	//VectorFliedd data3 (d1d, n, 1);
	switch(id){

		case 1: xCell = &d1d[0];
			n_xCell = n;
			//printf("xCell ");
			//for (int i=0; i<10; i++){
			//	printf("%f ", xCell[i]);
			//}
			break;

		case 2: yCell = &d1d[0];
			n_yCell = n;
			break;

		case 3: zCell = &d1d[0];
			n_zCell = n;
			break;

		case 4: 
			velocityX = &d1d[0];
			n_velocityX = n;
			//	printf("adapter2 velocityX ");
			//	for (int i=0;i<110;i++){
			//		printf("%f ", velocityX[i]);
			//	}
			break;

		case 5: velocityY = &d1d[0];
			n_velocityY = n;
			break;

		case 6: velocityZ = &d1d[0];
			n_velocityZ = n;
			break;

		case 7: zTop = &d1d[0];
			n_zTop = n;
			break;

		case 8: xVertex = &d1d[0];
			n_xVertex = n;
			//printf("xVertex[0] %f %f %f\n",xVertex[0], xVertex[1], xVertex[2] );
			break;

		case 9: yVertex = &d1d[0];
			n_yVertex = n;
			break;

		case 10: zVertex = &d1d[0];
			 n_zVertex = n;
			 break;

		default:
			 break;     
	}

	// if (num_arrived_fields == 2){
	//   decaf2->put(container);
	//	num_arrived_fields = 0;
	//   }
}

void func3(){
	printf("in func3\n");
	//decaf_put_int(1);
}

//void run (Workflow &workflow){

// MPI_Init(NULL, NULL);
// MPI_Comm C_communicator = MPI_Comm_f2c (communicator);
// printf(" before decaf obj construction\n");
// decaf2 = new Decaf(MPI_COMM_WORLD, workflow);
// printf("after decaf obj construction\n");
//if (decaf2->my_node("prod"))
//     prod();
// if (decaf2->my_node("con"))
//     con();

// cleanup
//delete decaf;


//}

void decaf_put_int(int n){

	SimpleFieldi data(n);
	pConstructData container;
	container->appendData("var", data,
			DECAF_NOFLAG, DECAF_PRIVATE,
			DECAF_SPLIT_KEEP_VALUE, DECAF_MERGE_ADD_VALUE);

	// send the data on all outbound dataflows
	// in this example there is only one outbound dataflow, but in general there could be more
	//decaf2->put(container, "out");

}


void finish_decaf(){
	if (decaf2->my_node("prod") || decaf2->my_node("con")){
		int pid=decaf2->world->rank();
		printf("pid %d terminating now \n", decaf2->world->rank());
		decaf2->terminate();
	}

	MPI_Finalize();
}



void f_MPI_Comm_c2f(MPI_Fint &comm) {

	comm  = decaf2->prod_comm_handle();

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



void prod(int n, double x[]){

	if (decaf2->my_node("prod")){
		printf("in producer\n");

		// the data in this example is just the timestep; add it to a container
		SimpleFieldi data(n);
		pConstructData container;
		container->appendData("var", data,
				DECAF_NOFLAG, DECAF_PRIVATE,
				DECAF_SPLIT_KEEP_VALUE, DECAF_MERGE_ADD_VALUE);

		// send the data on all outbound dataflows
		// in this example there is only one outbound dataflow, but in general there could be more
		decaf2->put(container, "out");
	}

}

void producer_terminate(){ 
	printf ("%d\n", decaf2->my_node("prod")); 
	if (decaf2->my_node("prod")){
		int pid=decaf2->world->rank();
		printf("producer %d terminating now \n", decaf2->world->rank());
		decaf2->terminate();
		printf("producer %d terminated \n", pid);
		MPI_Finalize();
	}else{
		printf("not in producer");
	}     

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
