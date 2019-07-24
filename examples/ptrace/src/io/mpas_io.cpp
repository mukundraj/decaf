#include "mpas_io.h"
#include "misc.h"

mpas_io::mpas_io(){

		velocityXv.resize(2);
		velocityYv.resize(2);
		velocityZv.resize(2);
		zTop.resize(2);
		zMid.resize(2);
		vertVelocityTop.resize(2);


		// TODO: replace hard code values, also by local sizes
		nCells = 7234;
		nEdges =  22736;
		nVertices = 15459;
		nVertLevels = 100;
		maxEdges = 6;
		nVertLevelsP1 = 101;


}

void mpas_io::update_data(int data_id, int frame_no, std::vector<int> &data_int, std::vector<double> &data_dbl){

	switch(data_id){

				case 0: //fprintf(stderr, "Recv xCell %d,\n", data_id);
						xCell = std::move(data_dbl);
						break;

				case 1: //fprintf(stderr, "Recv yCell %d,\n", data_id);
						yCell = std::move(data_dbl);
						break;

				case 2:	
						//printf(stderr, "Recv zCell %d,\n", data_id);
						zCell = std::move(data_dbl);
						break;

				case 3:	//fprintf(stderr, "Recv xVertex %d,\n", data_id);
						xVertex = std::move(data_dbl);
						break;
						
				case 4:	//fprintf(stderr, "Recv yVertex %d,\n", data_id);
						yVertex = std::move(data_dbl);
						break;

				case 5:	//fprintf(stderr, "Recv zVertex %d,\n", data_id);
						zVertex = std::move(data_dbl);
						break;

				case 6:	//fprintf(stderr, "Recv indexToVertexID %d,\n", data_id);
						indexToVertexID = std::move(data_int);
						break;

				case 7:	//fprintf(stderr, "Recv indexToCellID %d,\n", data_id);
						indexToCellID = std::move(data_int);

						for (int i=0; i<indexToCellID.size(); i++) { 
								cellIndex[indexToCellID[i]] = i;
						}
						break;

				case 8:	//fprintf(stderr, "Recv verticesOnEdge %d,\n", data_id);
						verticesOnEdge = std::move(data_int);
						break;

				case 9:	//fprintf(stderr, "Recv cellsOnVertex %d,\n", data_id);
						cellsOnVertex = std::move(data_int);
						break;

				case 10:{	//fprintf(stderr, "Recv verticesOnCell %d,\n", data_id);
						verticesOnCell = std::move(data_int);
						// size_t max=0;
						// 	for (size_t i=0; i<verticesOnCell.size(); i++){
						// 		if (max<verticesOnCell[i])
						// 			max = verticesOnCell[i];
						// 	}
							// dprint("verticesOnCell %ld", verticesOnCell.size());
						}
						break;
						

				case 11:	//fprintf(stderr, "Recv velocityXv %d, fno %d\n", data_id, frame_no);
						velocityXv[frame_no%2] = std::move(data_dbl);
						break;

				case 12:	//fprintf(stderr, "Recv velocityYv %d, fno %d\n", data_id, frame_no);
						velocityYv[frame_no%2] = std::move(data_dbl);
						break;

				case 13:	//fprintf(stderr, "Recv velocityZv %d,\n", data_id);
						velocityZv[frame_no%2] = std::move(data_dbl);
						break;

				case 14://	fprintf(stderr, "Recv nEdgesOnCell %d,\n", data_id);
						nEdgesOnCell = std::move(data_int);
						break;

				case 15://	fprintf(stderr, "Recv maxLevelCell %d,\n", data_id);
						maxLevelCell = std::move(data_int);
						break;

				case 16:	//fprintf(stderr, "Recv boundaryVertex %d,\n", data_id);
						boundaryVertex = std::move(data_int);
						break;

				case 17:	//fprintf(stderr, "Recv vertVelocityTop %d,\n", data_id);
						vertVelocityTop[frame_no%2] = std::move(data_dbl);
						break;

				case 18:{//fprintf(stderr, "Recv cellsOnCell %d,\n", data_id);
							cellsOnCell = std::move(data_int);
						} 
						
						break;

				case 19:	//fprintf(stderr, "Recv zTop %d,\n", data_id);
						zTop[frame_no%2] = std::move(data_dbl);
						break;

				case 20:	//fprintf(stderr, "Recv zMid %d,\n", data_id);
						zMid[frame_no%2] = std::move(data_dbl);
						break;

				default:fprintf(stderr, "Exiting! %d,\n", data_id); 
						exit(0);
						break;
				}

}