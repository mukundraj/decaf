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
#include "pblock.hpp"
#include "block_io.h"

#include <diy/mpi.hpp>
#include <diy/master.hpp>
#include <diy/assigner.hpp>
#include <diy/serialization.hpp>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "nabo/nabo.h"
#include "interpolators.h"

#include "misc.h"

//#include "advect.h"


// TODO: remove global variables
Eigen::MatrixXd *M_local;
Eigen::MatrixXd C_local;

void parallel_streamlines(mpaso &mpas_g, mpaso &mpas_c, int round, PBlock &b, int max_steps, const diy::Master::ProxyWithLink&   cp, const diy::Assigner& assigner, diy::mpi::communicator &world );

int frame_number = 0, global_gid; //FIX_ME: remove global variable

// computes vertex values for ghost vertices that also lie on native cells.
void compute_ghost_vertex_values(mpaso &mpas_c){

	int nCellsAll = mpas_c.xyzCell.size()/3;
	dprint("nCells nCellsAll %ld %d", mpas_c.nCells, nCellsAll);
	C_local.resize(3, nCellsAll);
	for (int i=0;i<nCellsAll;i++){
		C_local(0,i) = mpas_c.xyzCell[i*3];
		C_local(1,i) = mpas_c.xyzCell[i*3+1];
		C_local(2,i) = mpas_c.xyzCell[i*3+2];
	}

	mpas_c.nns_cells = Nabo::NNSearchD::createKDTreeLinearHeap(C_local);


	// TODO: 1) deal with pentagonal cells 2) Avoid repeated computations of velocity and top levels of vertices
	int nVertLevels = mpas_c.nVertLevels;

	mpas_c.velocityXv.resize(mpas_c.vertexIndex.size()*mpas_c.nVertLevels);
	mpas_c.velocityYv.resize(mpas_c.vertexIndex.size()*mpas_c.nVertLevels);
	mpas_c.velocityZv.resize(mpas_c.vertexIndex.size()*mpas_c.nVertLevels);
	mpas_c.zTopVertex.resize(mpas_c.vertexIndex.size()*mpas_c.nVertLevels);


	
	for (int i=0;i<mpas_c.nCells;i++){



		for(int j=0;j<6;j++){
			int nvgid = mpas_c.verticesOnCell[i*6+j]; // neighbor vertex global id
			// if(mpas_c.block_vert_ids.find(nvgid)==mpas_c.block_vert_ids.end() ){

				// get local vertex id from global id
			int nvlid = mpas_c.vertexIndex[nvgid];

			int c0 = mpas_c.cellIndex[mpas_c.cellsOnVertex[nvlid*3]], c1 = mpas_c.cellIndex[mpas_c.cellsOnVertex[nvlid*3+1]], c2 = mpas_c.cellIndex[mpas_c.cellsOnVertex[nvlid*3+2]];

				// if (mpas_c.indexToCellID[i]==4463){
				// 	dprint("4463 found in %d, nei nvgid vlid,  c0g c1g c2g, c0 c1 c2,  %d %d, %d %d %d, %d %d %d", mpas_c.gid, nvgid, nvlid, mpas_c.cellsOnVertex[nvlid*3], mpas_c.cellsOnVertex[nvlid*3+1], mpas_c.cellsOnVertex[nvlid*3+2], c0, c1, c2);
				// }



			double X[3][3] = {
				{mpas_c.xyzCell[c0*3], mpas_c.xyzCell[c0*3+1], mpas_c.xyzCell[c0*3+2]},
				{mpas_c.xyzCell[c1*3], mpas_c.xyzCell[c1*3+1], mpas_c.xyzCell[c1*3+2]},
				{mpas_c.xyzCell[c2*3], mpas_c.xyzCell[c2*3+1], mpas_c.xyzCell[c2*3+2]}
			};

			double P[3] = {mpas_c.xVertex[nvlid], mpas_c.yVertex[nvlid], mpas_c.zVertex[nvlid]};

			double lambda[3];
			barycentric_point2triangle(X[0], X[1], X[2], P, lambda);


			

			
			// if (0){
					// dprint("neighbor cell zero for %d", nvgid);

			for (int curVertLevel=0;curVertLevel<mpas_c.nVertLevels; curVertLevel++){

				mpas_c.zTopVertex[nvlid*mpas_c.nVertLevels + curVertLevel] = lambda[0] * mpas_c.zTop[c0*nVertLevels+curVertLevel] + lambda[1] * mpas_c.zTop[c1*nVertLevels+curVertLevel] + lambda[2] * mpas_c.zTop[c2*nVertLevels+curVertLevel];

				if (mpas_c.cellsOnVertex[nvlid*3]==0 || mpas_c.cellsOnVertex[nvlid*3+1] ==0 || mpas_c.cellsOnVertex[nvlid*3+2] == 0){
					mpas_c.velocityXv[nvlid*mpas_c.nVertLevels+curVertLevel] = 0;
					mpas_c.velocityYv[nvlid*mpas_c.nVertLevels+curVertLevel] = 0;
					mpas_c.velocityZv[nvlid*mpas_c.nVertLevels+curVertLevel] = 0;
					
				}else{

					for (int curVertLevel=0;curVertLevel<mpas_c.nVertLevels; curVertLevel++){



						mpas_c.velocityXv[nvlid*mpas_c.nVertLevels+curVertLevel] = lambda[0] * mpas_c.velocityX[c0*mpas_c.nVertLevels+curVertLevel] + lambda[1] * mpas_c.velocityX[c1*mpas_c.nVertLevels+curVertLevel] + lambda[2] * mpas_c.velocityX[c2*mpas_c.nVertLevels+curVertLevel];
						mpas_c.velocityYv[nvlid*mpas_c.nVertLevels+curVertLevel] = lambda[0] * mpas_c.velocityY[c0*mpas_c.nVertLevels+curVertLevel] + lambda[1] * mpas_c.velocityY[c1*mpas_c.nVertLevels+curVertLevel] + lambda[2] * mpas_c.velocityY[c2*mpas_c.nVertLevels+curVertLevel];
						mpas_c.velocityZv[nvlid*mpas_c.nVertLevels+curVertLevel] = lambda[0] * mpas_c.velocityZ[c0*mpas_c.nVertLevels+curVertLevel] + lambda[1] * mpas_c.velocityZ[c1*mpas_c.nVertLevels+curVertLevel] + lambda[2] * mpas_c.velocityZ[c2*mpas_c.nVertLevels+curVertLevel];
						// mpas_c.zTopVertex[nvlid*mpas_c.nVertLevels + curVertLevel] = lambda[0] * mpas_c.zTop[c0*nVertLevels+curVertLevel] + lambda[1] * mpas_c.zTop[c1*nVertLevels+curVertLevel] + lambda[2] * mpas_c.zTop[c2*nVertLevels+curVertLevel];

					}


				}
				
			}
			
			// }

		}


	}



}

void populate_block( PBlock* b,
	int gid,
	std::vector<double> &data_bar, mpaso *mpasobj, diy::mpi::communicator &world, mpaso &mpas_c ){

	int frame_no = (int) data_bar[14];
	int nVertices = (int) data_bar[15];
	int nCells = (int) data_bar[16];
	int nVertLevels = (int) data_bar[17];
	// int nCells_global;
	std::vector<int> nCell_counts, nVertices_counts; 
	std::vector<std::vector<double>> xCells_global, yCells_global, zCells_global, 
	xVertex_global, yVertex_global, zVertex_global; 

	int read_order[14] = {4, 5, 6, 11, 0, 1, 2, 3, 8, 9,10,7, 12, 13};
	int field_offsets[14];
	int field_sizes[14];
	int bsp = 0; // block start position

	field_offsets[read_order[0]] = bsp+20;
	field_sizes[read_order[0]] = data_bar[bsp+read_order[0]];
	for (int j=1; j<14; j++){
		field_offsets[read_order[j]] = field_offsets[read_order[j-1]]+field_sizes[read_order[j-1]];
		field_sizes[read_order[j]] = data_bar[bsp+read_order[j]];
	}

	fprintf(stderr, "inside populate %d %d %d %d %d\n", gid, frame_no, nVertices, nCells, nVertices);
	if (frame_no==1){
		b->indexToCellID.resize(nCells);  //0  
		b->xCell.resize(nCells);  //1
		b->yCell.resize(nCells); //2
		b->zCell.resize(nCells); //3

		b->velocityX.resize(nVertLevels*nCells); //4
		b->velocityY.resize(nVertLevels*nCells); //5
		b->velocityZ.resize(nVertLevels*nCells); //6

		b->indexToVertexID.resize(nVertices);	//7
		b->xVertex.resize(nVertices); //8
		b->yVertex.resize(nVertices); //9
		b->zVertex.resize(nVertices); //10
		b->zTop.resize(nVertLevels*nCells); //11
		b->cellsOnVertex.resize(nVertices*3); //12
		b->verticesOnCell.resize(nCells*6); //13

		// for (int i =0;i<13; i++)
		//   fprintf(stderr, "field offsets %d %d %d\n", i, field_offsets[i], field_sizes[i]);
		// populate frame 1 exclusive members 
		int k = 0; // xCell index offset var
		for (int j=field_offsets[0]; j<field_offsets[0]+field_sizes[0]; j++){
			b->xCell[k] = data_bar[field_offsets[1]+k]; // -1 on LHS to deal with fortran to C offset
			b->yCell[k] = data_bar[field_offsets[2]+k];
			b->zCell[k] = data_bar[field_offsets[3]+k];
			b->indexToCellID[k] = (int)data_bar[field_offsets[0]+k];
			
			// //std::cout<<" "<<int(data_bar[j]);
			k++;
		}

		k = 0;
		for (int j=field_offsets[7];j<field_offsets[7]+field_sizes[7];j++){
			b->indexToVertexID[k] = (int)data_bar[field_offsets[7]+k];
			b->xVertex[k] = data_bar[field_offsets[8]+k];
			b->yVertex[k] = data_bar[field_offsets[9]+k];
			b->zVertex[k] = data_bar[field_offsets[10]+k];

			k++;
		}

		int l=0;
		for (int j=field_offsets[7];j<field_offsets[7]+field_sizes[7];j++){
			// int m = data_bar[j] - 1;

			for (int k=0; k<3; k++){
				b->cellsOnVertex[l*3+k] = data_bar[field_offsets[12]+3*l+k];
			}
			l++;	
		}

		fprintf(stderr, "verticesOnCell %d %d %d\n",field_offsets[13], field_sizes[13], (int)data_bar.size());
		l=0;
		for (int j=field_offsets[0];j<field_offsets[0]+field_sizes[0];j++){
			for (int k=0; k<6; k++){
				b->verticesOnCell[l*6+k] = data_bar[field_offsets[13]+6*l+k];
			}
			l++;	
		}




		// MPI_allgather for getting number of cells
		diy::mpi::all_gather(world, nCells, nCell_counts );
		int nCells_global = 0;
		for (auto& n : nCell_counts)
			nCells_global += n;

		diy::mpi::all_gather(world, b->xCell, xCells_global);
		diy::mpi::all_gather(world, b->yCell, yCells_global);
		diy::mpi::all_gather(world, b->zCell, zCells_global);
		mpasobj->xCells.reserve(nCells_global);	
		mpasobj->yCells.reserve(nCells_global);	
		mpasobj->zCells.reserve(nCells_global);	
		for (size_t i=0;i<xCells_global.size();i++){
			mpasobj->xCells.insert(mpasobj->xCells.end(), xCells_global[i].begin(), xCells_global[i].end());

			mpasobj->yCells.insert(mpasobj->yCells.end(), yCells_global[i].begin(), yCells_global[i].end());

			mpasobj->zCells.insert(mpasobj->zCells.end(), zCells_global[i].begin(), zCells_global[i].end());
		}
		mpasobj->nCells = nCells_global;

		fprintf(stderr, "nCell_counts %ld %d %d\n", mpasobj->zCells.size(), nCells_global, nCell_counts[1]);

		// MPI_allgather for getting number of vertices
		diy::mpi::all_gather(world, nVertices, nVertices_counts);
		int nVertices_global = 0;
		for (auto &n : nVertices_counts)
			nVertices_global += n;

		diy::mpi::all_gather(world, b->xVertex, xVertex_global);
		diy::mpi::all_gather(world, b->yVertex, yVertex_global);
		diy::mpi::all_gather(world, b->zVertex, zVertex_global);
		mpasobj->xVertex.reserve(nVertices_global);	
		mpasobj->yVertex.reserve(nVertices_global);	
		mpasobj->zVertex.reserve(nVertices_global);	

		for (size_t i=0;i<xVertex_global.size();i++){
			mpasobj->xVertex.insert(mpasobj->xVertex.end(), xVertex_global[i].begin(), xVertex_global[i].end());

			mpasobj->yVertex.insert(mpasobj->yVertex.end(), yVertex_global[i].begin(), yVertex_global[i].end());

			mpasobj->zVertex.insert(mpasobj->zVertex.end(), zVertex_global[i].begin(), zVertex_global[i].end());
		}
		mpasobj->nVertices = nVertices_global;



		// Construct global kd trees for vertices and cells

		// Eigen::MatrixXd C(3, mpasobj->nCells);
		// for (int i=0;i<mpasobj->nCells;i++){
		// 	C(0,i) = mpasobj->xCells[i];
		// 	C(1,i) = mpasobj->yCells[i];
		// 	C(2,i) = mpasobj->zCells[i];
		// }

		// mpasobj->nns_cells = Nabo::NNSearchD::createKDTreeLinearHeap(C);

		// Eigen::MatrixXd M(3, mpasobj->nVertices);
		// for (int i=0;i<mpasobj->nVertices;i++){
		// 	M(0,i) = mpasobj->xVertex[i];
		// 	M(1,i) = mpasobj->yVertex[i];
		// 	M(2,i) = mpasobj->zVertex[i];
		// }


		// mpasobj->nns = Nabo::NNSearchD::createKDTreeLinearHeap(M);

		// Construct local kd trees for vertices and cells

		




	}
	// fprintf(stderr, "indexToCell %f %f\n", data_bar[field_offsets[0]+1], data_bar[field_offsets[0]+0] );

	fprintf(stderr, "xCell %f %f %d\n", b->xCell[0], b->xCell[1], b->indexToCellID[0] );

	

	// populate the recurring members	
	int l=0;
	for (int j=field_offsets[0]; j<field_offsets[0]+field_sizes[0];j++){
		// int m = data_bar[j]-1;
		// fprintf(stderr, " %d", m);
		for (int k=0;k<nVertLevels;k++){
			b->velocityX[l*nVertLevels+k] = data_bar[field_offsets[4]+nVertLevels*l+k];
			b->velocityY[l*nVertLevels+k] = data_bar[field_offsets[5]+nVertLevels*l+k];
			b->velocityZ[l*nVertLevels+k] = data_bar[field_offsets[6]+nVertLevels*l+k];
			b->zTop[l*nVertLevels+k] = data_bar[field_offsets[11]+nVertLevels*l+k];
			b->indexToCellID[l] = data_bar[j];
		}
		l++;
	}




	mpas_c.indexToCellID = b->indexToCellID;
	mpas_c.xCells = b->xCell;
	mpas_c.yCells = b->yCell;
	mpas_c.zCells = b->zCell;
	mpas_c.velocityX = b->velocityX;
	mpas_c.velocityY = b->velocityY;
	mpas_c.velocityZ = b->velocityZ;
	mpas_c.indexToVertexID = b->indexToVertexID;
	mpas_c.xVertex = b->xVertex;
	mpas_c.yVertex = b->yVertex;
	mpas_c.zVertex = b->zVertex;
	mpas_c.zTop = b->zTop;
	mpas_c.cellsOnVertex = b->cellsOnVertex;
	mpas_c.verticesOnCell = b->verticesOnCell;

	mpas_c.nCells = mpas_c.xCells.size();
	mpas_c.nVertices = mpas_c.xVertex.size();

	mpas_c.gid = b->gid;




	// C_local.resize(3, nCells);
	// for (int i=0;i<nCells;i++){
	// 	C_local(0,i) = mpas_c.xCells[i];
	// 	C_local(1,i) = mpas_c.yCells[i];
	// 	C_local(2,i) = mpas_c.zCells[i];
	// }

	// mpas_c.nns_cells = Nabo::NNSearchD::createKDTreeLinearHeap(C_local);

	M_local = new Eigen::MatrixXd(3, nVertices);
	for (int i=0;i<nVertices;i++){
		(*M_local)(0,i) = b->xVertex[i];
		(*M_local)(1,i) = b->yVertex[i];
		(*M_local)(2,i) = b->zVertex[i];
	}


	mpas_c.nns = Nabo::NNSearchD::createKDTreeLinearHeap(*M_local);

	
	if (frame_no==1 and b->gid==0){
		
		int lcell_0idx = 0;
		fprintf(stderr, "velocityX at %d and gid 0\n", lcell_0idx);
		for (int g=0;g<10;g++){
			fprintf(stderr, "%f ", mpas_c.velocityX[lcell_0idx*nVertLevels+g]);
		}
		fprintf(stderr, "\n");


		// fprintf(stderr, "verticesOnCell in main \n");
		// for (int i=0;i<10;i++)
		// 	fprintf(stderr, "%d\n", mpas_c.verticesOnCell[i]);
	}
	
	// int max_indexToVertexID=0;
	// fprintf(stderr, "indexToVertexID\n");
	// if (frame_no==1){
	// 	for (int g=0;g<nVertices;g++){
	// 		// fprintf(stderr, "%d ", mpas_c.indexToVertexID[g]);
	// 		if (max_indexToVertexID<mpas_c.indexToVertexID[g]){
	// 			max_indexToVertexID = mpas_c.indexToVertexID[g];
	// 		}

	// 	}
	// }
	// fprintf(stderr, "max_indexToVertexID %d\n", max_indexToVertexID);

	mpas_c.nVertLevels = 100; // FIX_ME: Remove the hard coding

	// for (int i=0; i<mpas_c.nCells; i++) { // shoudn't i vary from 1 to nCells here? then cellIndex[indexToCellID[i-1]] = i-1
	// 	mpas_c.cellIndex[mpas_c.indexToCellID[i]] = i;
	// 	// fprintf(stderr, "%d, %d\n", i, indexToCellID[i]);
	// 	//	std::cout<<i<<" "<<indexToCellID[i]<<" "<<cellIndex[i]<<"\n";
	// }



	// fprintf(stderr, "indexToCellID cellIndex\n");
	// if (frame_no==1){
	// 	for (int g=0;g<10;g++){
	// 		fprintf(stderr, "%d %d \n", mpas_c.indexToCellID[g], mpas_c.cellIndex[g]);
	// 	}
	// }
	// fprintf(stderr, "\n");

	// mapping from global to local vertex id
	for (int i=0; i<mpas_c.nVertices; i++) {
		mpas_c.vertexIndex[mpas_c.indexToVertexID[i]] = i;
		mpas_c.block_vert_ids.insert(mpas_c.indexToVertexID[i]);
		// fprintf(stderr, "%d, %d\n", i, indexToVertexID[i]);
	}


	

	// fprintf(stderr, "max_cell_on_vert %d\n", max_cell_on_vert);

	// if (b->gid==0){

	// 	for (int i=0;i<10;i++){
	// 		// fprintf(stderr, " xCell %d %f %f %f\n", i, mpas_c.xCells[i], mpas_c.yCells[i], mpas_c.zCells[i]);
	// 		fprintf(stderr, " xVertex %d %f %f %f\n", mpas_c.vertexIndex[i], mpas_c.xVertex[i], mpas_c.yVertex[i], mpas_c.zVertex[i]);
	// 	}
	// }

	// derive velocity on cell verticies 
	mpas_c.velocityXv.resize(mpas_c.nVertices*mpas_c.nVertLevels);
	mpas_c.velocityYv.resize(mpas_c.nVertices*mpas_c.nVertLevels);
	mpas_c.velocityZv.resize(mpas_c.nVertices*mpas_c.nVertLevels);
	mpas_c.zTopVertex.resize(mpas_c.nVertices*mpas_c.nVertLevels);
	mpas_c.zTopVertexNorm.resize(mpas_c.nVertices*mpas_c.nVertLevels);
	mpas_c.xyzCell.resize(mpas_c.nCells*3);
	for (int i=0; i<mpas_c.nCells; i++) mpas_c.xyzCell[i*3] = mpas_c.xCells[i];
		for (int i=0; i<mpas_c.nCells; i++) mpas_c.xyzCell[i*3+1] = mpas_c.yCells[i];
			for (int i=0; i<mpas_c.nCells; i++) mpas_c.xyzCell[i*3+2] = mpas_c.zCells[i];
				int max_cell_on_vert = 0;

	/* // commenting out vertex velocity and height computation as redone  in compute_ghost_vertex_values
	for (int i=0; i<mpas_c.nVertices; i++) {
		//    if (cellsOnVertex[i*3] == 0 || cellsOnVertex[i*3+1] == 0 || cellsOnVertex[i*3+2] == 0) continue; // on boundary
		// if (mpas_c.cellsOnVertex[i*3]>max_cell_on_vert)
		// 	max_cell_on_vert = mpas_c.cellsOnVertex[i*3];
		// if (mpas_c.cellsOnVertex[i*3+1]>max_cell_on_vert)
		// 	max_cell_on_vert = mpas_c.cellsOnVertex[i*3+1];
		// if (mpas_c.cellsOnVertex[i*3+2]>max_cell_on_vert)
		// 	max_cell_on_vert = mpas_c.cellsOnVertex[i*3+2];

				int c0 = mpas_c.cellIndex[mpas_c.cellsOnVertex[i*3]], c1 = mpas_c.cellIndex[mpas_c.cellsOnVertex[i*3+1]], c2 = mpas_c.cellIndex[mpas_c.cellsOnVertex[i*3+2]];



				double X[3][3] = {
					{mpas_c.xyzCell[c0*3], mpas_c.xyzCell[c0*3+1], mpas_c.xyzCell[c0*3+2]},
					{mpas_c.xyzCell[c1*3], mpas_c.xyzCell[c1*3+1], mpas_c.xyzCell[c1*3+2]},
					{mpas_c.xyzCell[c2*3], mpas_c.xyzCell[c2*3+1], mpas_c.xyzCell[c2*3+2]}
				};

				// if (i==2){
				// 	fprintf(stderr, "cellsOnVertex for i=2 %d %d %d %d %d %d %f %f %f %f %f %f  %f %f %f\n", mpas_c.cellsOnVertex[i*3], mpas_c.cellsOnVertex[i*3+1], mpas_c.cellsOnVertex[i*3+2], c0, c1, c2, mpas_c.xyzCell[c0*3], mpas_c.xyzCell[c0*3+1], mpas_c.xyzCell[c0*3+2], 
				// 		mpas_c.xyzCell[c1*3], mpas_c.xyzCell[c1*3+1], mpas_c.xyzCell[c1*3+2],
				// 		mpas_c.xyzCell[c2*3], mpas_c.xyzCell[c2*3+1], mpas_c.xyzCell[c2*3+2]
				// 		);
				// }

				// TODO: Check if the following lines can be removed since also being done in compute ghost vertex values

				double P[3] = {mpas_c.xVertex[i], mpas_c.yVertex[i], mpas_c.zVertex[i]};

				double lambda[3];
				barycentric_point2triangle(X[0], X[1], X[2], P, lambda);
		//[i * nVertLevels + curVertLevel];
				for (int curVertLevel=0;curVertLevel<mpas_c.nVertLevels; curVertLevel++){

					mpas_c.velocityXv[i*mpas_c.nVertLevels+curVertLevel] = lambda[0] * mpas_c.velocityX[c0*mpas_c.nVertLevels+curVertLevel] + lambda[1] * mpas_c.velocityX[c1*mpas_c.nVertLevels+curVertLevel] + lambda[2] * mpas_c.velocityX[c2*mpas_c.nVertLevels+curVertLevel];
					mpas_c.velocityYv[i*mpas_c.nVertLevels+curVertLevel] = lambda[0] * mpas_c.velocityY[c0*mpas_c.nVertLevels+curVertLevel] + lambda[1] * mpas_c.velocityY[c1*mpas_c.nVertLevels+curVertLevel] + lambda[2] * mpas_c.velocityY[c2*mpas_c.nVertLevels+curVertLevel];
					mpas_c.velocityZv[i*mpas_c.nVertLevels+curVertLevel] = lambda[0] * mpas_c.velocityZ[c0*mpas_c.nVertLevels+curVertLevel] + lambda[1] * mpas_c.velocityZ[c1*mpas_c.nVertLevels+curVertLevel] + lambda[2] * mpas_c.velocityZ[c2*mpas_c.nVertLevels+curVertLevel];

					mpas_c.zTopVertex[mpas_c.nVertLevels*i + curVertLevel] = lambda[0] * mpas_c.zTop[c0*nVertLevels+curVertLevel] + lambda[1] * mpas_c.zTop[c1*nVertLevels+curVertLevel] + lambda[2] * mpas_c.zTop[c2*nVertLevels+curVertLevel];

					// if (i==2 && curVertLevel==50){
					// fprintf(stderr, "cur velocity level 50 at vertex for i=2 %f %f %f\n", 
					// 	mpas_c.velocityXv[i*mpas_c.nVertLevels+curVertLevel], 
					// 	mpas_c.velocityYv[i*mpas_c.nVertLevels+curVertLevel],
					// 	mpas_c.velocityZv[i*mpas_c.nVertLevels+curVertLevel]
					// 	);

					// 	fprintf(stderr, "velocityX l=50 i=2 %f\n", mpas_c.velocityX[c0*mpas_c.nVertLevels+curVertLevel]);
					// }


				}
	}
	*/


		}



		void enq_halo_info(mpaso &mpas_g, mpaso &mpas_c, int cur_bid, int nprocs, 
			PBlock*                              b,
			const diy::Master::ProxyWithLink&   cp,
			const diy::Assigner& assigner){

	std::vector<std::set<int>> ghost_cell_ids; // 

	// std::vector<ghost_req> ghost_cell_reqs; // 
	
	

	ghost_cell_ids.resize(nprocs);
	
	// ghost_cell_reqs.resize(nprocs);




	// populate vector of required cell ids

	for (int i=0; i<mpas_c.nCells; i++){
		int cur_cell_1gid = mpas_c.indexToCellID[i];

		for (int j=0; j<mpas_g.cell_g_neighbors[cur_cell_1gid].size(); j++){
			
			// check if the current neighbor is in same partition/block
			int nbr_cell_1gid = mpas_g.cell_g_neighbors[cur_cell_1gid][j];

			int cell_bid = mpas_g.cgid_to_bid[nbr_cell_1gid];

			if (cell_bid != cur_bid){
				ghost_cell_ids[cell_bid].insert(nbr_cell_1gid);
				

			}

		}



	}
	

	// remote queue
	for (int i=0;i<ghost_cell_ids.size();i++){

		if (ghost_cell_ids[i].size()>0){
			// fprintf(stderr, "\n");
			// for ( auto it = ghost_cell_ids[i].begin(); it != ghost_cell_ids[i].end(); it++ )
			// 	fprintf(stderr, "%d, ",*it);
			// fprintf(stderr, "\n");

			std::vector<int> ghost_cell_ids_vec;
			for ( auto it = ghost_cell_ids[i].begin(); it != ghost_cell_ids[i].end(); it++ )
				ghost_cell_ids_vec.push_back(*it);

			// ghost_cell_reqs[i].requester = cur_bid;
			fprintf(stderr, "ghost cell id size %ld\n", ghost_cell_ids[i].size());
			ghost_cell_ids_vec.push_back(cur_bid); // putting the current block id as last element in list

			int dest_gid            = i;
			int dest_proc           = assigner.rank(dest_gid);
			diy::BlockID dest_block = {dest_gid, dest_proc};
			cp.enqueue(dest_block, ghost_cell_ids_vec);
			// cp.enqueue(dest_block, 1 );
			fprintf(stderr, "enque info %d %d %d\n", i, dest_gid, dest_proc);

		}
		
	}


}

void process_halo_req(PBlock* b, const diy::Master::ProxyWithLink& cp, const diy::Assigner& assigner, mpaso mpas_c){

	vector<int> in;
	

	cp.incoming(in);
	
	fprintf(stderr, "in size %ld\n", in.size());
	for (size_t i = 0; i < in.size(); i++)
		if (cp.incoming(in[i]).size())
		{	
			std::vector<int> ghost_cell_ids;
			
			cp.dequeue(in[i],ghost_cell_ids);
			fprintf(stderr, "dequeued %ld\n", ghost_cell_ids.size());
			int dest_gid = ghost_cell_ids.back();
			ghost_cell_ids.pop_back();

			// for ( auto it = ghost_cell_ids.begin(); it != ghost_cell_ids.end(); it++ )
			// 	fprintf(stderr, "%d ", *it);

			// TODO Deal with local variables inside loops
			// gather cell info
			std::vector<ghost_point> ghost_info;
			for (int j=0;j<ghost_cell_ids.size();j++){
				ghost_point gp;



				int cellIdxLocal = mpas_c.cellIndex[ghost_cell_ids[j]];

				// include this cell info in ghost point
				gp.cgid = ghost_cell_ids[j];
				int lgid = mpas_c.cellIndex[gp.cgid];
				gp.cx = mpas_c.xyzCell[3*lgid];
				gp.cy = mpas_c.xyzCell[3*lgid+1];
				gp.cz = mpas_c.xyzCell[3*lgid+2];



				// gp.vel_cx.insert(gp.vel_cx.end(), &mpas_c.velocityX[lgid*mpas_c.nVertLevels], &mpas_c.velocityX[(1+lgid)*mpas_c.nVertLevels]);
				// gp.vel_cy.insert(gp.vel_cy.end(), &mpas_c.velocityY[lgid*mpas_c.nVertLevels], &mpas_c.velocityY[(1+lgid)*mpas_c.nVertLevels]);
				// gp.vel_cz.insert(gp.vel_cz.end(), &mpas_c.velocityZ[lgid*mpas_c.nVertLevels], &mpas_c.velocityZ[(1+lgid)*mpas_c.nVertLevels]);
				for (int k=0;k<mpas_c.nVertLevels;k++){
					gp.vel_cx[k] = mpas_c.velocityX[lgid*mpas_c.nVertLevels + k];
					gp.vel_cy[k] = mpas_c.velocityY[lgid*mpas_c.nVertLevels + k];
					gp.vel_cz[k] = mpas_c.velocityZ[lgid*mpas_c.nVertLevels + k];
					gp.zTop[k] = mpas_c.zTop[lgid*mpas_c.nVertLevels + k];
				}



				int ngvid=-1; // neighbor global vertex id
				// get neighboring vertices
				for (int k=0;k<6;k++){
					// TODO: investigate crash on setting following to -1
					gp.vert_gids[k] = -1; // initially assume vertex not present
					
					
					ngvid = mpas_c.verticesOnCell[lgid*6+k];


					// if (b->gid==1){
					// 	fprintf(stderr, "spacking vertex %d %d\n", gp.cgid, ngvid);

					// 	} 
					

					// if vertex present in partition. include this vertex's info in ghost point.
					if (mpas_c.vertexIndex.find(ngvid)!=mpas_c.vertexIndex.end()){
						


						
						int nlvid = mpas_c.vertexIndex[ngvid]; // get local vid

						gp.vert_gids[k] = ngvid; // store global vertex id
						gp.vert_x[k] = mpas_c.xVertex[nlvid];
						gp.vert_y[k] = mpas_c.yVertex[nlvid];
						gp.vert_z[k] = mpas_c.zVertex[nlvid];

						// populate vertex neighbor cell-ids
						gp.cellsOnVertex[k*3] = mpas_c.cellsOnVertex[nlvid*3];
						gp.cellsOnVertex[k*3+1] = mpas_c.cellsOnVertex[nlvid*3+1];
						gp.cellsOnVertex[k*3+2] = mpas_c.cellsOnVertex[nlvid*3+2];

						// if (ngvid==14425){
						// 		fprintf(stderr, "found inend1 14425, %d %d %f\n", k, ngvid, gp.vert_x[k]);
						// 		}

					}

				}

				// push gp to ghost_info
				ghost_info.push_back(gp);


			}

			
			




			// queue cell info for exchange
			
			// if (cp.gid()==0)
			// 	dest_gid=1;
			// else
			// 	dest_gid=0;
			int dest_proc           = assigner.rank(dest_gid);
			diy::BlockID dest_block = {dest_gid, dest_proc};
			cp.enqueue(dest_block, ghost_info);


		}

		dprint("end of process process_halo_req in gid %d", b->gid);

	}

	void deq_halo_info(PBlock* b, const diy::Master::ProxyWithLink& cp, mpaso &mpas_c){
		vector<int> in;


		cp.incoming(in);

		fprintf(stderr, "in size %ld\n", in.size());
		for (size_t i = 0; i < in.size(); i++){
			std::vector<ghost_point> ghost_info;
			if (cp.incoming(in[i]).size())
			{	
				cp.dequeue(in[i],ghost_info);
			// fprintf(stderr, "dequeued again %ld %d %f %f\n", ghost_info.size(), ghost_info[0].cgid, ghost_info[0].cx, ghost_info[0].vel_cx[1]);

			// append xCell, yCell, zCell, associated vertices pos + indices, create ghost index map
				for (int j=0;j<ghost_info.size(); j++ ){

				// process the cell
					mpas_c.xyzCell.push_back(ghost_info[j].cx);
					mpas_c.xyzCell.push_back(ghost_info[j].cy);
					mpas_c.xyzCell.push_back(ghost_info[j].cz);
					std::copy(ghost_info[j].vel_cx, ghost_info[j].vel_cx+mpas_c.nVertLevels, std::back_inserter(mpas_c.velocityX));
					std::copy(ghost_info[j].vel_cy, ghost_info[j].vel_cy+mpas_c.nVertLevels, std::back_inserter(mpas_c.velocityY));
					std::copy(ghost_info[j].vel_cz, ghost_info[j].vel_cz+mpas_c.nVertLevels, std::back_inserter(mpas_c.velocityZ));
					std::copy(ghost_info[j].zTop, ghost_info[j].zTop+mpas_c.nVertLevels, std::back_inserter(mpas_c.zTop));

				// update indexToCellID, cellIndex
					mpas_c.indexToCellID.push_back(ghost_info[j].cgid);
					mpas_c.cellIndex[ghost_info[j].cgid] = mpas_c.indexToCellID.size()-1;

					if (ghost_info[j].cgid == 4462 || ghost_info[j].cgid == -1){
						dprint("found cgid %d in ghost in block %d", ghost_info[j].cgid, b->gid);
						dprint("velX at 70 %f", mpas_c.velocityX[mpas_c.nVertLevels*ghost_info[j].cgid+70]);
						dprint("velY at 70 %f", mpas_c.velocityY[mpas_c.nVertLevels*ghost_info[j].cgid+70]);
						dprint("velZ at 70 %f", mpas_c.velocityZ[mpas_c.nVertLevels*ghost_info[j].cgid+70]);
					}	

				// update ghost cell ids into std::set

					

					for (int k=0;k<6;k++){
						if (ghost_info[j].vert_gids[k]>-1){
						// process associated vertex


							

						// vertexIndex
							int vgid = ghost_info[j].vert_gids[k]; 

							if (mpas_c.vertexIndex.count(vgid) == 0){

								
								
								mpas_c.vertexIndex[vgid] = (int) mpas_c.xVertex.size();

								if (ghost_info[j].vert_gids[k]==14425){
									dprint("found inend 14425, %d %d %f %f %f %d %d %d", mpas_c.vertexIndex[vgid], k, ghost_info[j].vert_x[k], ghost_info[j].vert_y[k], ghost_info[j].vert_z[k], ghost_info[j].cellsOnVertex[k*3], ghost_info[j].cellsOnVertex[k*3+1], ghost_info[j].cellsOnVertex[k*3]+2);
								}

								// // vertex positions
								mpas_c.xVertex.push_back(ghost_info[j].vert_x[k]);
								mpas_c.yVertex.push_back(ghost_info[j].vert_y[k]);
								mpas_c.zVertex.push_back(ghost_info[j].vert_z[k]);

								// // cellOnVertex?*not needed for pure ghost vertices. needed for shared ghost. keep for all.
								mpas_c.cellsOnVertex.push_back(ghost_info[j].cellsOnVertex[k*3]);
								mpas_c.cellsOnVertex.push_back(ghost_info[j].cellsOnVertex[k*3+1]);
								mpas_c.cellsOnVertex.push_back(ghost_info[j].cellsOnVertex[k*3+2]);

							} 


						}
					}

				}


			}
		}

	// compute velocities at shared vertices after halo exchange. Only native and shared.
	}


// enqueue remote data
// there is still a link, but you can send to any BlockID = (gid, proc)
	void remote_enq(
		PBlock*                              b,
		const diy::Master::ProxyWithLink&   cp,
		const diy::Assigner&                assigner)
	{
	// as a test, send my gid to block 1 gids away (in a ring), which is outside the link
	// (the link has only adjacent block gids)
		int my_gid              = cp.gid();
		int dest_gid            = (my_gid + 1) % assigner.nblocks();
		int dest_proc           = assigner.rank(dest_gid);
		diy::BlockID dest_block = {dest_gid, dest_proc};

	// b->values[0] = my_gid;

	// cp.enqueue(dest_block, my_gid);

		std::vector<EndPt> endptvec;
		endptvec.resize(3);

		cp.enqueue(dest_block, endptvec);
	}

// dequeue remote data
// there is still a link, but exchange(remote = true) exchanged messages from any block
	void remote_deq(PBlock* b, const diy::Master::ProxyWithLink& cp)
	{
		std::vector<int> incoming_gids;
		cp.incoming(incoming_gids);
		for (size_t i = 0; i < incoming_gids.size(); i++)
			if (cp.incoming(incoming_gids[i]).size())
			{
				vector<EndPt> recvd_data;
				cp.dequeue(incoming_gids[i], recvd_data);
				fmt::print(stderr, "Remote dequeue: gid {} received value {} from gid {}\n",
					cp.gid(), recvd_data[0].pid, incoming_gids[i]);
				fprintf(stderr, "Received vector: %d\n", 3);
			}
		}






// consumer
		void con(Decaf* decaf)
		{      	
			int max_rounds = 3;
	// int max_steps = 2000;

			int max_steps = 400;

			mpaso mpas_c, mpas_g;
			vector< pConstructData > in_data;


			diy::mpi::communicator world(decaf->con_comm_handle());
			fprintf(stderr, "diy world size: %d", world.size() );

			diy::FileStorage storage("./DIY.XXXXXX");
			int                       nblocks    = world.size();
			int                       threads    = 1;
			int                       mem_blocks = -1;

			diy::Master               master(world,
				threads,
				mem_blocks,
				&create_block,
				&destroy_block,
				&storage,
				&save_block,
				&load_block);


			diy::RoundRobinAssigner assigner(world.size(), nblocks);

	// this example creates a linear chain of blocks with links between adjacent blocks
	std::vector<int> gids;                     // global ids of local blocks
	assigner.local_gids(world.rank(), gids);   // get the gids of local blocks

	PBlock *b;
	int pgid;
	block_io bio;
	for (size_t i = 0; i < gids.size(); ++i)   // for the local blocks in this processor
	{
		int gid = gids[i];

		diy::Link*   link = new diy::Link;   // link is this block's neighborhood
		//        diy::BlockID neighbor;               // one neighbor in the neighborhood
		//        if (gid < nblocks - 1)               // RH neighbors for all but the last block in the global domain
		//        {
		//            neighbor.gid  = gid + 1;                     // gid of the neighbor block
		//            neighbor.proc = assigner.rank(neighbor.gid); // process of the neighbor block
		//            link->add_neighbor(neighbor);                // add the neighbor block to the link
		//        }
		//        if (gid > 0)                         // LH neighbors for all but the first block in the global domain
		//        {
		//            neighbor.gid  = gid - 1;
		//            neighbor.proc = assigner.rank(neighbor.gid);
		//            link->add_neighbor(neighbor);
		//        }

		// keep pointer to current block. use it for populating / updating only current block later
		b = new PBlock;
		pgid = gid;
		global_gid = gid;
		b->gid = gid; // TODO: remove other gid vars e.g. global_gid
		master.add(gid, b, link);    // add the current local block to the master
	}


	// if (b->gid==0){
	// 	int i = 0;

	// 	char hostname[256];
	// 	gethostname(hostname, sizeof(hostname));
	// 	printf("PID %d on %s ready for attach\n", getpid(), hostname);
	// 	fflush(stdout);
	// 	while (0 == i)
	// 		sleep(5);
	// }


	std::string ip_file = "graph.topology";
	mpas_g.generate_domain_decomposition_graph(ip_file, world.size());	
	mpas_c.cgid_to_bid = mpas_g.cgid_to_bid; // TODO: compute this val directly in mpas_c/object b
	mpas_g.read_cell_g_neighbors();
	// (remote) exchange some data outside the links
	//	master.foreach([&](PBlock* b, const diy::Master::ProxyWithLink& cp)
	//		{ 	
	//			remote_enq(b, cp, assigner); }); // seed/deque particles, update vector fields, trace the epoch, and enqueue unfinished particles


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
				// printf("got data_bar in con: %f %f %f %f %f %f %f %f %f %f\n", data_bar[13], data_bar[14], data_bar[15], data_bar[16], data_bar[17], data_bar[18], data_bar[1452848+0], data_bar[1452848+4], data_bar[1452848+5], data_bar[1452848+6]);
				// fir_cell_idx = data_bar[0];
				// //printf("fir cell val %f\n", data_bar[int(12+data_bar[4]+data_bar[5]+data_bar[6]+data_bar[7])]);
				// fprintf(stderr, "databar size :%ld\n", data_bar.size());
				// unsigned int init_t=0, fin_t = 40000000, h = 200000;
				// std::vector<double> seeds_xyz(1000*3);
				// std::string op_file = std::to_string(ctr)+".vtp";
				//streamlines slines(mpas1, data_bar, op_file , init_t, fin_t, h, seeds_xyz);

				// use assigner to
				// master.foreach([&](PBlock* b, const diy::Master::ProxyWithLink& cp)
				//                 {
				//                 populate_block(b, cp, assigner, data_bar ); }); // seed/deque particles, update vector fields, trace the epoch, and enqueue unfinished particles
				int frame_no = (int) data_bar[14];

				if (frame_no==1){

					// streamlines slines;
					mpas_c.compute_cellIndex(b->gid, world.size());
					// TODO  premature stopping if work completed before max_rounds 
					populate_block(b, pgid, data_bar, &mpas_g, world, mpas_c);

					// identify and enque halo cells
					master.foreach([&](PBlock* b, const diy::Master::ProxyWithLink& cp)
						{ enq_halo_info(mpas_g, mpas_c, b->gid, world.size(), 
							b, cp, assigner); });
					
					//rexchange
					bool remote = true;
					master.exchange(remote);
					// process the halo request
					master.foreach([&](PBlock* b, const diy::Master::ProxyWithLink& cp){
						process_halo_req(b, cp, assigner, mpas_c);
					});
					// rexchange
					master.exchange(remote);

					//deq_halo_info
					master.foreach([&](PBlock* b, const diy::Master::ProxyWithLink& cp){
						deq_halo_info(b, cp, mpas_c);
					});

					// if (frame_no==1 && b->gid==0){

					// 	int cid = mpas_c.cellIndex[4462];

					// 	fprintf(stderr, "verticesOnCell 4462\n");

					// 	for (int u=2;u<3;u++){
					// 		int lvid = mpas_c.vertexIndex[mpas_c.verticesOnCell[cid*6+u]];
					// 		fprintf(stderr, "%d %d %f %f %f, ", mpas_c.verticesOnCell[cid*6+u], lvid, mpas_c.xVertex[lvid], mpas_c.yVertex[lvid], mpas_c.zVertex[lvid]);
					// 	}
					// 	dprint("\nprinting vertex local ids (vlids)");
					// 	for (int u=0;u<6;u++){
					// 		fprintf(stderr, "%d, ", mpas_c.vertexIndex[mpas_c.verticesOnCell[cid*6+u]]);
					// 	}
					// 	fprintf(stderr, "\n");


					// }

					compute_ghost_vertex_values(mpas_c);

					int stop = (max_rounds ? max_rounds : 1);
					int incr = (max_rounds ? 1 : 0);
					for (int round = 0; round < stop; round += incr)
					{

						master.foreach([&](PBlock* b, const diy::Master::ProxyWithLink& cp)
						{

							parallel_streamlines(mpas_g, mpas_c, round, *b, max_steps, cp,  assigner, world);
						});
						bool remote = true;
						master.exchange(remote);	
					}

					//bio.write_cell_centers(pgid, world, b->xCell); //only testing
					bio.write_particle_traces(pgid, world, *b, max_steps, mpas_c);
				}	

			}else{ 
				fprintf(stderr, "null ptr data_bar\n");
			}


			// foreach load data into block


			// another loop for block synchronous tracing, just keep one traceblock functioni
			// followed by exchange	


			// bool remote = true;
			// master.exchange(remote);
			// master.foreach(&remote_deq);	// particle  tracing here
		}
		// printf("consumer sum = %d\n", sum);
		ctr++;
	}


	// terminate the task (mandatory) by sending a quit message to the rest of the workflow
	fprintf(stderr, "consumer terminating\n");
	decaf->terminate();

}

int main(){



	dprint("starting mpas consumer");



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
	delete M_local;
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

//returns global cell id or -1 if out of vertically out of global
bool get_curpos_vel_sl(mpaso &mpas_g, mpaso &mpas_c, double *X, Eigen::Vector3d &c_vel, bool &in_global_domain){

	int K = 6;
	double depth = std::sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]) - mpas_g.radius;

	// if (mpas_c.gid==1){
	// 	dprint("depth 1: %f", depth);
	// }
	// if (mpas_c.gid==0){
	// 	dprint("depth 0: %f", depth);
	// }

	// get nearest vertex neighbor ids using mpas_c
	Eigen::VectorXi nearest_idx(K);
	Eigen::VectorXd dists2(K);

	Eigen::VectorXd q(3);
	q<<X[0], X[1], X[2];
	// q << -5682112.802109 , -619878.622472 , 2814587.637284;

	

	// mpas_c.nns->knn(q, nearest_idx, dists2, K, 0, Nabo::NNSearchF::SORT_RESULTS);
	

	// get nearest cell neighbor ids using mpas_c
	Eigen::VectorXi nearest_cell_idx(1);
	Eigen::VectorXd dists2_cell(1);
	mpas_c.nns_cells->knn(q, nearest_cell_idx,dists2_cell, 1,0, Nabo::NNSearchF::SORT_RESULTS);

	int clid = nearest_cell_idx[0]; // cell local id
	for (int i=0;i<K;i++){
		int gvid = mpas_c.verticesOnCell[clid*6+i];
		nearest_idx[i] = mpas_c.vertexIndex[gvid];
	}
	// dprint("clid %d %d %f %f %f", clid, mpas_c.indexToCellID[clid], q[0], q[1], q[2]);
	if (clid>mpas_c.nCells-1){
		return false; // out of partition
	}

	

	// if (mpas_c.gid == 1){
	// 	dprint("nearest_cell_idx for 1: %d %d, %d %d %d %d %d %d", clid, mpas_c.indexToCellID[clid], 
	// 		mpas_c.vertexIndex[mpas_c.verticesOnCell[clid*6+0]],
	// 		mpas_c.vertexIndex[mpas_c.verticesOnCell[clid*6+1]],
	// 		mpas_c.vertexIndex[mpas_c.verticesOnCell[clid*6+2]],
	// 		mpas_c.vertexIndex[mpas_c.verticesOnCell[clid*6+3]],
	// 		mpas_c.vertexIndex[mpas_c.verticesOnCell[clid*6+4]],
	// 		mpas_c.vertexIndex[mpas_c.verticesOnCell[clid*6+5]]);

	// }

	// if (mpas_c.gid == 0){
	// 	dprint("X %f %f %f", X[0], X[1], X[2]);
	// 	dprint("nearest_cell_idx for gid0 %d %f %f %f %d %f %f %f", nearest_cell_idx[0], 
	// 		mpas_c.xCells[nearest_cell_idx[0]], mpas_c.yCells[nearest_cell_idx[0]], mpas_c.zCells[nearest_cell_idx[0]], mpas_c.cellIndex[4462], mpas_c.xCells[mpas_c.cellIndex[4462]], mpas_c.yCells[mpas_c.cellIndex[4462]], mpas_c.zCells[mpas_c.cellIndex[4462]]);

	// }



	// if (mpas_c.gid == 0){
	// 	for (int i=0;i<K;i++){
	// 	// nearest_idx_array[i] = nearest_idx[i];
	// 	// fprintf(stderr, "nearest %d %f %f %f %f\n", nearest_idx_array[i], mpas_c.xVertex[i], 
	// 	// 	mpas_c.yVertex[i], mpas_c.zVertex[i], dists2[i]);

	// 	// nearest_cell_idx -> local. get global cell idx. global vert neighbors, local vert neighbors.
	// 		int idx = nearest_cell_idx[0];
	// 		int gcid = mpas_c.indexToCellID[idx];
	// 		int gvid = mpas_c.verticesOnCell[idx*6+i];
	// 		int lvid = mpas_c.vertexIndex[gvid];

	// 		dprint("nverts %d %d %d %d %f %f %f", idx, gcid, gvid, lvid, mpas_c.xVertex[lvid], mpas_c.yVertex[lvid], mpas_c.zVertex[lvid]);

	// 	}
	// }


	// interpolate vertically and get positions and velocity (in values)
	std::vector<double> values;
	values.resize(6*nearest_idx.size());


	double top_boun = mpas_c.zTop[mpas_c.nVertLevels*nearest_cell_idx[0]+0];

	double botom_boun = mpas_c.zTop[mpas_c.nVertLevels*nearest_cell_idx[0]+99];


	if (depth> top_boun || depth< botom_boun){
		// terminate flow line
		in_global_domain = false;
		dprint("vertical exit printing values: %f %f %f", top_boun, botom_boun, depth);
		return false; 
		
	}else{ 
		interpolate_vertically(mpas_c.nVertLevels, mpas_c.zTopVertex, nearest_idx, values, depth, mpas_c.xVertex, mpas_c.yVertex, mpas_c.zVertex, mpas_c.velocityXv, mpas_c.velocityYv, mpas_c.velocityZv, mpas_c.gid);

		// interpolate horizontally
		interpolate_horizontally(X[0], X[1], X[2], values, c_vel);
		// if (!interpolate_horizontally(X[0], X[1], X[2], values, c_vel))
		// 	return false;
	}

	if (std::isnan(c_vel[0]) || std::isnan(c_vel[1]) || std::isnan(c_vel[2])){
		c_vel[0] = values[0];
		c_vel[1] = values[1];
		c_vel[2] = values[2];
	}
	int cgid = mpas_c.indexToCellID[clid];
	if (mpas_c.gid == 0 && cgid == 6){
		dprint("nearest_cell_idx for 0 clid cgid: %d %d, %d %d %d %d %d %d, %f %f %f", clid, cgid, nearest_idx[0], nearest_idx[1], nearest_idx[2], nearest_idx[3], nearest_idx[4], nearest_idx[5], c_vel[0], c_vel[1], c_vel[2] );
	}

	// dprint("interpolated velocity %d: %f %f %f %f %f %f", mpas_c.gid, X[0], X[1], X[2], c_vel[0], c_vel[1], c_vel[2]);

	// if (global_gid==1 ){
		// fprintf(stderr,"check vel %f %f %f %f %f %f\n", mpas_c.velocityXv[5539],
		//                                          mpas_c.velocityXv[909],
		//                                           mpas_c.velocityXv[5565],
		//                                            mpas_c.velocityXv[1899],
		//                                             mpas_c.velocityXv[1880],
		//                                              mpas_c.velocityXv[7875]);

	


		//cout<<"inside  get_vel "<<nearest_idx[0]<<" "<<nearest_idx[1]<<" "<<nearest_idx[2]<<" "<<nearest_idx[3]<<" "<<nearest_idx[4]<<" "<<nearest_idx[5]<<" "<<c_vel<<" "<<global_gid<<" "<<mpas_c.velocityXv[5539]<<" "<<mpas_c.velocityXv[909]<<" "<<mpas_c.velocityXv[5565]<<" "<<mpas_c.velocityXv[1899]<<" "<<mpas_c.velocityXv[1880]<<" "<<mpas_c.velocityXv[7875]<<"\n";
	// }
	return true; // True implies curpoint within partition and within global domain vertically

}

Eigen::MatrixXd project_on_sphere(double *X, double *Y, double *P){

    // const double PI = std::acos(-1.0); // TODO: Put PI in a common namespace

	double xy[3] = {Y[0]-X[0], Y[1]-X[1], Y[2]-X[2]};
	double hx = std::sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
	double d = std::sqrt((xy[0])*(xy[0]) + (xy[1])*(xy[1]) + (xy[2])*(xy[2]));
	double hy = std::sqrt(Y[0]*Y[0] + Y[1]*Y[1] + Y[2]*Y[2]);

	double theta = acos((-d*d + hx*hx + hy*hy)/(2*hx*hy));

	double phi = d/hx;

//    fprintf(stderr, "hts %f %f %f %f %f\n", hx, d, hy, theta*180/PI, phi*180/PI);

    // get axis
	double a[3];
	cross_product(X, xy, a);
	double a_m = std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	a[0] /= a_m;
	a[1] /= a_m;
	a[2] /= a_m;

//    fprintf(stderr, "axis %f %f %f\n",a[0], a[1], a[2]);

	auto angle = phi;
	auto sinA = std::sin(angle / 2);
	auto cosA = std::cos(angle / 2);

	Eigen::Quaterniond q;
	q.x() = a[0] * sinA;
	q.y() = a[1] * sinA;
	q.z() = a[2] * sinA;
	q.w() = cosA;
	Eigen::MatrixXd R = q.normalized().toRotationMatrix();
	Eigen::Vector3d x(X[0], X[1], X[2]);
	Eigen::Vector3d p = R*x;

	P[0] = p[0];
	P[1] = p[1];
	P[2] = p[2]; 

	return R;

}

double mag(double *x){
	return std::sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}


bool rk1_mpas(mpaso &mpas_g,
	mpaso &mpas_c,
	double *X,
	double h,
	double *Y,
	bool &in_global_domain, 
	PBlock &b)
{

	
	Eigen::Vector3d c_vel;

	if (get_curpos_vel_sl(mpas_g, mpas_c, X, c_vel, in_global_domain)){	

		// advect
		Y[0] = X[0] + h*c_vel[0];	
		Y[1] = X[1] + h*c_vel[1];
		Y[2] = X[2] + h*c_vel[2];
		// dprint("magnitudes %f %f, %f %f %f, %f %f %f", mag(X), mag(Y), X[0], X[1], X[2], Y[0], Y[1], Y[2]);

		Eigen::MatrixXd R = project_on_sphere(X, Y, Y);

		// if (b.gid==0){
		// 	dprint("magnitudes %f %f, %f %f %f, %f, %d", mag(X), mag(Y), X[0], X[1], X[2], R.determinant(), R.transpose()*R==R.transpose()*R);
		// }

		return true;

	}

	return false;
	// else{
	// 	in_global_domain = false;
	// 	return false; // needs to be true for continued advection
	// }

}

void parallel_streamlines(mpaso &mpas_g, mpaso &mpas_c, int round, PBlock &b, int max_steps, const diy::Master::ProxyWithLink&   cp, const diy::Assigner& assigner, diy::mpi::communicator &world){

	// fprintf(stderr, "inside parallel %f\n", mpas_c.xCell_p[0]);
	vector<EndPt> particles;
	map<diy::BlockID, vector<EndPt> > outgoing_endpts;

	double dc_t, dc_c, ratio, dc_b, dc_base; //dist from center of earth of top and current layer

	
	// if first iteration, then seed

	if (round==0){
		// int skipval = 5000;
		int skipval = 2500;
		// std::vector<int> seed_level_ids = {10, 30, 50, 70, 90};
		std::vector<int> seed_level_ids = { 10};
		int n_seed_verts = mpas_c.nVertices/skipval;
		int n_seed_levels = seed_level_ids.size();

		int n_particles = n_seed_verts * n_seed_levels;
		diy::mpi::all_gather(world, n_particles, b.global_nP);
		b.global_start_ids.resize(b.global_nP.size());
		// fprintf(stderr, "global np size %ld \n", b.global_nP.size());
		
		int n_particles_global = b.global_nP[0];
		for (int k = 1; k<b.global_nP.size(); k++){
			b.global_start_ids[k] = b.global_start_ids[k-1]+b.global_nP[k-1];
			// fprintf(stderr, "global_start_pos %d %d\n", k, b.global_start_pos[k]);
			n_particles_global += b.global_nP[k];
		}
		b.global_trace_sizes.resize(n_particles_global);
		dprint("n_seed_verts %d b.gid %d", n_seed_verts, b.gid);
		for (int i=0;i<n_seed_verts; i++){

			for (int j=0;j<n_seed_levels;j++){

				int vert_id = i*skipval;
				// if (b.gid==1){
				// 	dprint("gid 1 seed_level_ids[j] %d zTopVertex %f",seed_level_ids[j],mpas_c.zTopVertex[mpas_c.nVertLevels*vert_id+seed_level_ids[j]]);
				// }


				dc_c = mpas_c.radius + mpas_c.zTopVertex[mpas_c.nVertLevels*vert_id+seed_level_ids[j]];
				dc_t = mpas_c.radius + mpas_c.zTopVertex[mpas_c.nVertLevels*vert_id+0];
				dc_b = mpas_c.radius + mpas_c.zTopVertex[mpas_c.nVertLevels*vert_id+99];


				// fprintf(stderr, "ratio %f\n", ratio);
				// ratio = dc_c/dc_t;
				// ratio = dc_c/mpas_c.radius;

				// EndPt q;
				// q[0] = mpas_c.xVertex[vert_id];
				// q[1] = mpas_c.yVertex[vert_id];
				// q[2] = mpas_c.zVertex[vert_id];

				dc_base = std::sqrt(mpas_c.xVertex[vert_id]*mpas_c.xVertex[vert_id] + mpas_c.yVertex[vert_id]*mpas_c.yVertex[vert_id] + mpas_c.zVertex[vert_id]*mpas_c.zVertex[vert_id]);
				ratio = dc_c/dc_base;

				EndPt p;
				p.pid = mpas_c.init;
				p.sid = mpas_c.init;
				p.step = 0;
				p.gpid = b.global_start_ids[b.gid] + p.pid;
				p[0] = mpas_c.xVertex[vert_id] * ratio;
				p[1] = mpas_c.yVertex[vert_id] * ratio;
				p[2] = mpas_c.zVertex[vert_id] * ratio;



				// if (b.gid==0){
				// 	p[0] = -2089982.455199 * ratio;
				// 	p[1] = 3032836.244409 * ratio;
				// 	p[2] = -5198695.665453 * ratio;
				// 	dprint("seed %f %f %f", p[0], p[1], p[2]); 

				// }

				// if (mpas_c.gid==1){
				// 	dprint("depth during init gid 1: %f %f, dc_c, dc_t, dc_b ratio %f %f %f %f",std::sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]) - mpas_g.radius, std::sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]) - mpas_g.radius, dc_c, dc_t, dc_b, ratio);
				// }

				dprint("vlid gvid bid %d %d %d", vert_id, mpas_c.indexToVertexID[vert_id], b.gid);
				// int gvid = mpas_c.indexToVertexID[vert_id];
				// dprint("cellsOnVertex %d %d %d %d", gvid, mpas_c.cellsOnVertex[vert_id*3], mpas_c.cellsOnVertex[vert_id*3+1], mpas_c.cellsOnVertex[vert_id*3+2]);
				particles.push_back(p);
				mpas_c.init++;



			}
		}

		

	}

	
	// get incoming endpts 	
	vector<int> incoming_gids;
	cp.incoming(incoming_gids);
	for (size_t i = 0; i < incoming_gids.size(); i++)
	{
		if (cp.incoming(incoming_gids[i]).size())
		{
			vector<EndPt> incoming_endpts;
			cp.dequeue(incoming_gids[i], incoming_endpts);
			for (size_t j = 0; j < incoming_endpts.size(); j++) {
				incoming_endpts[j].sid++;
				particles.push_back(incoming_endpts[j]);
				dprint("received particle");
			}


		}

	}	


	// trace particles
	int nenq_particles = 0;

	for (int i = 0; i < particles.size(); i++)
	{
		Pt&     cur_p = particles[i].pt; // current end point
		Segment s(particles[i]);         // segment with one point p (initialized with EndPt)
		Pt      next_p;                  // coordinates of next end point
		bool finished = false, in_global_domain = true;
		// int cur_start = s.start_step;
		// int cur_segment_steps = 0;
		// trace particle as far as it will go
		
		while(rk1_mpas(mpas_g, mpas_c, cur_p.coords, 100000, next_p.coords, in_global_domain, b)){	
		// while(rk1_mpas(mpas_g, mpas_c, cur_p.coords, 10000, next_p.coords, in_global_domain)){	
			s.pts.push_back(next_p);
			cur_p = next_p;
			// if (s.pts.size()%100==0)
				// dprint(" %d %ld", b.gid, s.pts.size()  );

			// dprint("inloop %d %d %ld %d %d %d", b.gid, mpas_c.get_bid_for_pt(next_p.coords), s.pts.size(), s.pid, s.sid, s.start_step);
			// cur_segment_steps++;
			// if (s.pts.size() >= max_steps) // CHECK cur_step here instead of s.pts.size()?

			// TODO: is the nan check really needed (deal with nan earlier skip_val=3000)?
			if (s.start_step + s.pts.size() >= max_steps || std::isnan(next_p.coords[0]) || std::isnan(next_p.coords[1]) || std::isnan(next_p.coords[2])) 
			{
				dprint("Reached max steps");
				finished = true;
				break;

			} 


		}
		// dprint("clid %d %d %f %f %f", b.gid, mpas_c.get_bid_for_pt(next_p.coords), next_p.coords[0], next_p.coords[1], next_p.coords[2]);
		
		// fprintf(stderr, "pid size %d %ld\n", s.pid, b.trace_sizes.size());
		// b.trace_sizes[s.pid] = s.start_step + s.pts.size();
		b.global_trace_sizes[s.gpid] = s.start_step + s.pts.size();

		// fprintf(stderr, "pts.size %ld\n", s.pts.size());

		if (in_global_domain==false)
			finished = true;


		b.segments.push_back(s);	

		if (finished){
			mpas_c.done++;

			dprint("finished particle in %d", b.gid);
		}
		else{ 			// package unfinished endpoint for sending 
			dprint("enqueuing particle");
			EndPt out_pt(s); // TO_DO: Also update the start_time here
			int dest_gid = mpas_c.get_bid_for_pt(next_p.coords);		
			dprint("b.gid dest_gid %d %d", b.gid, dest_gid);
			int dest_proc = assigner.rank(dest_gid);

			diy::BlockID bid = {dest_gid, dest_proc};

			outgoing_endpts[bid].push_back(out_pt);


			nenq_particles++;

		}	

	}	
	

	// send off unfinished segments using rexchange

	for (map<diy::BlockID, vector<EndPt> >::const_iterator it =
		outgoing_endpts.begin(); it != outgoing_endpts.end(); it++)
		cp.enqueue(it->first, it->second);

	

}


void block_io::write_particle_traces(int gid, const diy::mpi::communicator &world, PBlock &b, int max_steps, mpaso &mpas_c){

	int ret, ncfile, nprocs, rank;
	char filename[256], buf[13] = "Hello World\n";
	int data[2];


	// int nP = b.trace_sizes.size(); // number of particles in current block
	// std::vector<std::vector<int>> trace_sizes_global;
	// concatenate trace_from all blocs in proc that processes block 0

	// diy::mpi::gather(world, b.trace_sizes, trace_sizes_global, 0);

	// if (rank == 0){
	//         diy::mpi::reduce(world, b.global_trace_sizes, b.global_trace_sizes,0, MPI_MAX);
	// }
	// else{
	//         diy::mpi::reduce(world, b.global_trace_sizes);
	// }
	//
	// // fprintf(stderr, "np %d\n", nP);
	//
	//

	MPI_Comm_rank(world, &rank);
	MPI_Comm_size(world, &nprocs);
	std::vector<int> trace_sizes_recv;
	trace_sizes_recv.resize(b.global_trace_sizes.size());
	// fprintf(stderr, "trace_sizes_global size %ld %d\n", trace_sizes_global.size(), rank);

	int trace_sizes_global[b.global_trace_sizes.size()];
	for (int i=0; i<b.global_trace_sizes.size();i++){
		trace_sizes_global[i] = b.global_trace_sizes[i];
	}	

	int global_max_segments = b.segments.size();

	if (gid==0){	
		MPI_Reduce(MPI_IN_PLACE, &trace_sizes_global, b.global_trace_sizes.size(), MPI_INT, MPI_MAX, 0, world);
		for (int i=0; i<b.global_trace_sizes.size();i++){
			dprint("trace_sizes_global i %d %d",i, trace_sizes_global[i]);
		}
		fprintf(stderr, " \n");
		MPI_Allreduce(MPI_IN_PLACE, &global_max_segments, 1, MPI_INT, MPI_MAX, world);
	}else{

		MPI_Reduce(&trace_sizes_global, &trace_sizes_global, b.global_trace_sizes.size(), MPI_INT, MPI_MAX, 0, world);
		MPI_Allreduce(MPI_IN_PLACE, &global_max_segments, 1, MPI_INT, MPI_MAX, world);
	}

	// dprint("global_max_segments %d gid %d", global_max_segments, b.gid );

	strcpy(filename, "particle_traces.nc");


	int dimid_s, dimid_p; // varid for steps and particles
	int ndims_pos = 2;
	int varid_xPos, varid_yPos, varid_zPos, varid_trace_sizes_global;
	int pos_dims[ndims_pos];

	ret = ncmpi_create(world, filename,
		NC_CLOBBER|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncfile);
	if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

	// define dimensions
	ret = ncmpi_def_dim(ncfile, "nStep", NC_UNLIMITED, &dimid_s);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	ret = ncmpi_def_dim(ncfile, "nParticle", b.global_trace_sizes.size(), &dimid_p);
	// ret = ncmpi_def_dim(ncfile, "nParticle", 4, &dimid_p); 
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	pos_dims[0] = dimid_s;
	pos_dims[1] = dimid_p;
	// define variables
	ret = ncmpi_def_var(ncfile, "xPos", NC_DOUBLE, ndims_pos, pos_dims, &varid_xPos);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	ret = ncmpi_def_var(ncfile, "yPos", NC_DOUBLE, ndims_pos, pos_dims, &varid_yPos);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	ret = ncmpi_def_var(ncfile, "zPos", NC_DOUBLE, ndims_pos, pos_dims, &varid_zPos);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	// MPI_Offset start[ndims_pos], count[ndims_pos];

	// end define mode
	ret = ncmpi_enddef(ncfile);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	// write data to file each segment as a separate call

	// TODO create and allocate buffer 
	// http://cucis.ece.northwestern.edu/projects/PnetCDF/doc/pnetcdf-c/ncmpi_005fput_005fvarn_005f_003ctype_003e.html
	MPI_Offset **starts, **counts;
	
	int num_reqs = b.segments.size();
	/* allocate starts and counts */


	starts    = (MPI_Offset**) malloc(num_reqs*       sizeof(MPI_Offset*));
	starts[0] = (MPI_Offset*)  calloc(num_reqs*ndims_pos, sizeof(MPI_Offset));
	for (int i=1; i<num_reqs; i++)
		starts[i] = starts[i-1] + ndims_pos;

	counts    = (MPI_Offset**) malloc(num_reqs*       sizeof(MPI_Offset*));
	counts[0] = (MPI_Offset*)  calloc(num_reqs*ndims_pos, sizeof(MPI_Offset));
	for (int i=1; i<num_reqs; i++)
		counts[i] = counts[i-1] + ndims_pos;
	
	// for (int i=0; i<num_reqs; i++){
	// 	starts[i][0] = b.segments[i].start_step; starts[i][1] = b.segments[i].gpid;
	// 	counts[i][0] = b.segments[i].pts.size(); counts[i][1] = 1;
	// 	dprint("segment gid %d gpid %d len %ld start_step %d", b.gid, b.segments[i].gpid, b.segments[i].pts.size(), b.segments[i].start_step);

	// 	// if (b.segments[i].gpid==1){
	// 	// 	dprint("here how %d %d %ld", b.gid, rank, b.segments.size());
	// 	// 	// for (int q=0;q<b.segments[i].pts.size();q++)
	// 	// 	// 	dprint("seg5 %d %d %f ", b.gid, q, b.segments[i].pts[q].coords[0]);
	// 	// }
	// }

	for (int i=0; i<num_reqs; i++){
		starts[i][0] = b.segments[i].start_step; starts[i][1] = b.segments[i].gpid;
		counts[i][0] = b.segments[i].pts.size(); counts[i][1] = 1;
		dprint("segment gid %d gpid %d len %ld start_step %d", b.gid, b.segments[i].gpid, b.segments[i].pts.size(), b.segments[i].start_step);

		
	}

       	/* allocate write buffer */
	int buf_len = 0;
	for (int i=0; i<num_reqs; i++) {
		MPI_Offset w_req_len=1;
		for (int j=0; j<ndims_pos; j++)
			w_req_len *= counts[i][j];
		buf_len += w_req_len;
	}


	dprint( "buflen num_reqs, gid %d %d %d", buf_len, num_reqs, b.gid);
	double *buffer = (double*) malloc(buf_len * sizeof(double));
	double *buffer_y = (double*) malloc(buf_len * sizeof(double));
	double *buffer_z = (double*) malloc(buf_len * sizeof(double));

	// for (int i=0; i<buf_len; i++) {
	//         buffer[i] = (double)rank;
	//         fprintf(stderr, "bufval %f\n", buffer[i]);
	// }
	
	// int bidx = 0;
	// for (int i=0; i<b.segments.size(); i++) {
	// 	for (int j=0; j<b.segments[i].pts.size();j++){
	// 		buffer[bidx] = b.segments[i].pts[j].coords[0];
	// 		buffer_y[bidx] = b.segments[i].pts[j].coords[1];
	// 		buffer_z[bidx] = b.segments[i].pts[j].coords[2];
	// 		bidx++;
	// 	}
	// }        

	int bidx = 0;
	for (int i=0; i<b.segments.size(); i++) {
		for (int j=0; j<b.segments[i].pts.size();j++){
			buffer[bidx] = b.segments[i].pts[j].coords[0];
			buffer_y[bidx] = b.segments[i].pts[j].coords[1];
			buffer_z[bidx] = b.segments[i].pts[j].coords[2];
			bidx++;
		}
	}        

	dprint("after allocation. %d %d", num_reqs, bidx);
	// dprint("after allocation. %d %f", bidx, buffer[13]); 	

	MPI_Offset start[2], count[2];
	
	// http://cucis.ece.northwestern.edu/projects/PnetCDF/
    // double *buffer_test = (double*) malloc(1 * sizeof(double));

	for (int i=0; i<global_max_segments; i++){

		
		dprint("segments size %ld %d %d", b.segments.size(), b.gid, i);
		if (i<b.segments.size()){
			double buffer_test[b.segments[i].pts.size()][1];

			for (int j=0;j<b.segments[i].pts.size();j++)
				buffer_test[j][0] = b.segments[i].pts[j].coords[0];

			start[0] = b.segments[i].start_step; start[1] = b.segments[i].gpid;
			count[0] = b.segments[i].pts.size(); count[1] = 1;
			ret = ncmpi_put_vara_double_all(ncfile, varid_xPos, start, count, &buffer_test[0][0]);
			if (ret != NC_NOERR) handle_error(ret, __LINE__);

			for (int j=0;j<b.segments[i].pts.size();j++)
				buffer_test[j][0] = b.segments[i].pts[j].coords[1];
			ret = ncmpi_put_vara_double_all(ncfile, varid_yPos, start, count, &buffer_test[0][0]);
			if (ret != NC_NOERR) handle_error(ret, __LINE__);

			for (int j=0;j<b.segments[i].pts.size();j++)
				buffer_test[j][0] = b.segments[i].pts[j].coords[2];
			ret = ncmpi_put_vara_double_all(ncfile, varid_zPos, start, count, &buffer_test[0][0]);
			if (ret != NC_NOERR) handle_error(ret, __LINE__);
			
		}else{

			double buffer_test[1][1];

			start[0] = 0; start[1] = 0;
			count[0] = 0; count[1] = 0;
			ret = ncmpi_put_vara_double_all(ncfile, varid_xPos, start, count, &buffer_test[0][0]);
			if (ret != NC_NOERR) handle_error(ret, __LINE__);
			ret = ncmpi_put_vara_double_all(ncfile, varid_yPos, start, count, &buffer_test[0][0]);
			if (ret != NC_NOERR) handle_error(ret, __LINE__);
			ret = ncmpi_put_vara_double_all(ncfile, varid_zPos, start, count, &buffer_test[0][0]);
			if (ret != NC_NOERR) handle_error(ret, __LINE__);
			
		}
		// dprint("written %d %d", i, b.gid);

	}


 //    double buffer_test[1][1];
 //    buffer_test[0][0] = 89;
 //    start[0] = 0; start[1] = rank;
 //    count[0] = 1; count[1] = 1;
	// ret = ncmpi_put_vara_double_all(ncfile, varid_xPos, start, count, &buffer_test[0][0]);
	// if (ret != NC_NOERR) handle_error(ret, __LINE__);

	// buffer_test[0][0] = 99;
 //    start[0] = rank; start[1] = 0;
 //    count[0] = 1; count[1] = 1;
	// ret = ncmpi_put_vara_double_all(ncfile, varid_xPos, start, count, &buffer_test[0][0]);
	// if (ret != NC_NOERR) handle_error(ret, __LINE__);

	// free(buffer_test);

	// ret = ncmpi_put_varn_double_all(ncfile, varid_xPos, num_reqs, starts, counts, buffer);
	// if (ret != NC_NOERR) handle_error(ret, __LINE__);
	// dprint("after first"); 

	// ret = ncmpi_put_varn_double_all(ncfile, varid_yPos, num_reqs ,starts, counts, buffer_y);
	// if (ret != NC_NOERR) handle_error(ret, __LINE__);

	// ret = ncmpi_put_varn_double_all(ncfile, varid_zPos, num_reqs ,starts, counts, buffer_z);
	// if (ret != NC_NOERR) handle_error(ret, __LINE__);

	for (int i=0;i<b.global_trace_sizes.size();i++){
		dprint("global %d %d", trace_sizes_global[i], b.gid);	
	}
	
	// close file
	ret = ncmpi_close(ncfile);
	if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

	dprint("after close"); 

	free(buffer);
	free(buffer_y);
	free(buffer_z);
	free(starts[0]);
	free(starts);
	free(counts[0]);
	free(counts);

	
	
	if (b.gid == 0){



		int dimid_nVertices, dimid_nVertexAllLayers, varid_xVertex, 
		varid_yVertex, varid_zVertex, 
		varid_velocityXv, varid_velocityYv, varid_velocityZv, varid_zTopVertex;
		size_t nVertices = mpas_c.xVertex.size();
		size_t nVertexAllLayers = 100 * nVertices;

		// create buffer arrays
		double *xVertex, *xVelocityXv;
		xVertex = new double [nVertices];
		xVelocityXv = new double [nVertexAllLayers];




		ret = ncmpi_open(MPI_COMM_SELF, filename,
			NC_WRITE|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncfile);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		
		ret = ncmpi_redef(ncfile); /* enter define mode */
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		ret = ncmpi_def_dim(ncfile, "nVertices", nVertices, &dimid_nVertices);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		
		ret = ncmpi_def_dim(ncfile, "nVertexAllLayers", nVertexAllLayers, &dimid_nVertexAllLayers);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		const size_t start_vertices[1] = {0}, size_vertices[1] = {nVertices};
		ret = ncmpi_def_var(ncfile, "xVertex", NC_DOUBLE, 1, &dimid_nVertices, &varid_xVertex);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);
		ret = ncmpi_def_var(ncfile, "yVertex", NC_DOUBLE, 1, &dimid_nVertices, &varid_yVertex);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);
		ret = ncmpi_def_var(ncfile, "zVertex", NC_DOUBLE, 1, &dimid_nVertices, &varid_zVertex);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);


		
		const size_t start_velocityV[1] = {0}, size_velocityV[1] = {nVertexAllLayers};
		ret = ncmpi_def_var(ncfile, "velocityVx", NC_DOUBLE, 1, &dimid_nVertexAllLayers, &varid_velocityXv);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		ret = ncmpi_def_var(ncfile, "velocityVy", NC_DOUBLE, 1, &dimid_nVertexAllLayers, &varid_velocityYv);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		ret = ncmpi_def_var(ncfile, "velocityVz", NC_DOUBLE, 1, &dimid_nVertexAllLayers, &varid_velocityZv);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		ret = ncmpi_def_var(ncfile, "zTopVertex", NC_DOUBLE, 1, &dimid_nVertexAllLayers, &varid_zTopVertex);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);
		
		//	Next define trace_sizes global
		ret = ncmpi_def_var(ncfile, "trace_sizes_global", NC_INT, 1, &dimid_p, &varid_trace_sizes_global);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		


		// end define mode
		ret = ncmpi_enddef(ncfile);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);
		
		

		MPI_Offset start[1], count[1];
		start[0]=0, count[0]=b.global_trace_sizes.size();
		ret = ncmpi_put_vara_int_all(ncfile, varid_trace_sizes_global, start, count, trace_sizes_global);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);	

		// xVertex copy to buffer arrays and write variable 
		for (size_t i=0;i<nVertices;i++){
			xVertex[i] = mpas_c.xVertex[i];
		}
		for (size_t i=0;i<nVertexAllLayers;i++){
			xVelocityXv[i] = mpas_c.velocityXv[i];
		}

		start[0]=0, count[0]=nVertices;
		ret = ncmpi_put_vara_double_all(ncfile, varid_xVertex, start, count, xVertex);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);
		start[0]=0, count[0]=nVertexAllLayers;
		ret = ncmpi_put_vara_double_all(ncfile, varid_velocityXv, start, count, xVelocityXv);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		// yVertex copy to buffer arrays and write variable
		for (size_t i=0;i<nVertices;i++){
			xVertex[i] = mpas_c.yVertex[i];
		}
		for (size_t i=0;i<nVertexAllLayers;i++){
			xVelocityXv[i] = mpas_c.velocityYv[i];
		}

		start[0]=0, count[0]=nVertices;
		ret = ncmpi_put_vara_double_all(ncfile, varid_yVertex, start, count, xVertex);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);
		start[0]=0, count[0]=nVertexAllLayers;
		ret = ncmpi_put_vara_double_all(ncfile, varid_velocityYv, start, count, xVelocityXv);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		// zVertex copy to buffer arrays and write variable
		for (size_t i=0;i<nVertices;i++){
			xVertex[i] = mpas_c.zVertex[i];
		}
		for (size_t i=0;i<nVertexAllLayers;i++){
			xVelocityXv[i] = mpas_c.velocityZv[i];
		}

		start[0]=0, count[0]=nVertices;
		ret = ncmpi_put_vara_double_all(ncfile, varid_zVertex, start, count, xVertex);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);
		start[0]=0, count[0]=nVertexAllLayers;
		ret = ncmpi_put_vara_double_all(ncfile, varid_velocityZv, start, count, xVelocityXv);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);


		for (size_t i=0;i<nVertexAllLayers;i++){
			xVelocityXv[i] = mpas_c.zTopVertex[i];
		}
		start[0]=0, count[0]=nVertexAllLayers;
		ret = ncmpi_put_vara_double_all(ncfile, varid_zTopVertex, start, count, xVelocityXv);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);
		



		// close file
		ret = ncmpi_close(ncfile);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		// free buffer arrays
		delete [] xVertex;
		delete [] xVelocityXv;


	}

	if (b.gid == 1){

		dprint("in process 1");

		int dimid_nVertices, dimid_nVertexAllLayers, varid_xVertex, 
		varid_yVertex, varid_zVertex, 
		varid_velocityXv, varid_velocityYv, varid_velocityZv, varid_zTopVertex;
		size_t nVertices = mpas_c.xVertex.size();
		size_t nVertexAllLayers = 100 * nVertices;

		// create buffer arrays
		double *xVertex, *xVelocityXv;
		xVertex = new double [nVertices];
		xVelocityXv = new double [nVertexAllLayers];


		strcpy(filename, "particle_traces2.nc");

		ret = ncmpi_create(MPI_COMM_SELF, filename,
			NC_WRITE|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncfile);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		
		// ret = ncmpi_redef(ncfile); /* enter define mode */
		// if (ret != NC_NOERR) handle_error(ret, __LINE__);

		ret = ncmpi_def_dim(ncfile, "nVertices", nVertices, &dimid_nVertices);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		
		ret = ncmpi_def_dim(ncfile, "nVertexAllLayers", nVertexAllLayers, &dimid_nVertexAllLayers);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		const size_t start_vertices[1] = {0}, size_vertices[1] = {nVertices};
		ret = ncmpi_def_var(ncfile, "xVertex", NC_DOUBLE, 1, &dimid_nVertices, &varid_xVertex);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);
		ret = ncmpi_def_var(ncfile, "yVertex", NC_DOUBLE, 1, &dimid_nVertices, &varid_yVertex);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);
		ret = ncmpi_def_var(ncfile, "zVertex", NC_DOUBLE, 1, &dimid_nVertices, &varid_zVertex);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);


		
		const size_t start_velocityV[1] = {0}, size_velocityV[1] = {nVertexAllLayers};
		ret = ncmpi_def_var(ncfile, "velocityVx", NC_DOUBLE, 1, &dimid_nVertexAllLayers, &varid_velocityXv);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		ret = ncmpi_def_var(ncfile, "velocityVy", NC_DOUBLE, 1, &dimid_nVertexAllLayers, &varid_velocityYv);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		ret = ncmpi_def_var(ncfile, "velocityVz", NC_DOUBLE, 1, &dimid_nVertexAllLayers, &varid_velocityZv);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		ret = ncmpi_def_var(ncfile, "zTopVertex", NC_DOUBLE, 1, &dimid_nVertexAllLayers, &varid_zTopVertex);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);
		
		// //	Next define trace_sizes global
		// ret = ncmpi_def_var(ncfile, "trace_sizes_global", NC_INT, 1, &dimid_p, &varid_trace_sizes_global);
		// if (ret != NC_NOERR) handle_error(ret, __LINE__);

		


		// end define mode
		ret = ncmpi_enddef(ncfile);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);
		
		

		MPI_Offset start[1], count[1];
		// start[0]=0, count[0]=b.global_trace_sizes.size();
		// ret = ncmpi_put_vara_int_all(ncfile, varid_trace_sizes_global, start, count, trace_sizes_global);
		// if (ret != NC_NOERR) handle_error(ret, __LINE__);	

		// xVertex copy to buffer arrays and write variable 
		for (size_t i=0;i<nVertices;i++){
			xVertex[i] = mpas_c.xVertex[i];
		}
		for (size_t i=0;i<nVertexAllLayers;i++){
			xVelocityXv[i] = mpas_c.velocityXv[i];
		}

		start[0]=0, count[0]=nVertices;
		ret = ncmpi_put_vara_double_all(ncfile, varid_xVertex, start, count, xVertex);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);
		start[0]=0, count[0]=nVertexAllLayers;
		ret = ncmpi_put_vara_double_all(ncfile, varid_velocityXv, start, count, xVelocityXv);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		// yVertex copy to buffer arrays and write variable
		for (size_t i=0;i<nVertices;i++){
			xVertex[i] = mpas_c.yVertex[i];
		}
		for (size_t i=0;i<nVertexAllLayers;i++){
			xVelocityXv[i] = mpas_c.velocityYv[i];
		}

		start[0]=0, count[0]=nVertices;
		ret = ncmpi_put_vara_double_all(ncfile, varid_yVertex, start, count, xVertex);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);
		start[0]=0, count[0]=nVertexAllLayers;
		ret = ncmpi_put_vara_double_all(ncfile, varid_velocityYv, start, count, xVelocityXv);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		// zVertex copy to buffer arrays and write variable
		for (size_t i=0;i<nVertices;i++){
			xVertex[i] = mpas_c.zVertex[i];
		}
		for (size_t i=0;i<nVertexAllLayers;i++){
			xVelocityXv[i] = mpas_c.velocityZv[i];
		}

		start[0]=0, count[0]=nVertices;
		ret = ncmpi_put_vara_double_all(ncfile, varid_zVertex, start, count, xVertex);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);
		start[0]=0, count[0]=nVertexAllLayers;
		ret = ncmpi_put_vara_double_all(ncfile, varid_velocityZv, start, count, xVelocityXv);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);


		for (size_t i=0;i<nVertexAllLayers;i++){
			xVelocityXv[i] = mpas_c.zTopVertex[i];
		}
		start[0]=0, count[0]=nVertexAllLayers;
		ret = ncmpi_put_vara_double_all(ncfile, varid_zTopVertex, start, count, xVelocityXv);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);
		



		// close file
		ret = ncmpi_close(ncfile);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		// free buffer arrays
		delete [] xVertex;
		delete [] xVelocityXv;


	}
}
