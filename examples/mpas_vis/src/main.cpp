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
#include "pblock.h"
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

void parallel_streamlines( int round, PBlock &b, int max_steps, const diy::Master::ProxyWithLink&   cp, const diy::Assigner& assigner, diy::mpi::communicator &world );

int frame_number = 0, global_gid; //FIX_ME: remove global variable

// computes vertex values for ghost vertices that also lie on native cells.
void compute_ghost_vertex_values(PBlock &b){

	int nCellsAll = b.xyzCell.size()/3;
	dprint("nCells nCellsAll %ld %d", b.nCells, nCellsAll);
	C_local.resize(3, nCellsAll);
	for (int i=0;i<nCellsAll;i++){
		C_local(0,i) = b.xyzCell[i*3];
		C_local(1,i) = b.xyzCell[i*3+1];
		C_local(2,i) = b.xyzCell[i*3+2];
	}

	b.nns_cells = Nabo::NNSearchD::createKDTreeLinearHeap(C_local);
	// mpas_c.nns_cells = Nabo::NNSearchD::createKDTreeLinearHeap(C_local);


	// TODO: 1) deal with pentagonal cells 2) Avoid repeated computations of velocity and top levels of vertices
	int nVertLevels = b.nVertLevels;

	// mpas_c.velocityXv.resize(mpas_c.vertexIndex.size()*mpas_c.nVertLevels);
	// mpas_c.velocityYv.resize(mpas_c.vertexIndex.size()*mpas_c.nVertLevels);
	// mpas_c.velocityZv.resize(mpas_c.vertexIndex.size()*mpas_c.nVertLevels);
	// mpas_c.zTopVertex.resize(mpas_c.vertexIndex.size()*mpas_c.nVertLevels);

	b.velocityXv.resize(b.vertexIndex.size()*b.nVertLevels);
	b.velocityYv.resize(b.vertexIndex.size()*b.nVertLevels);
	b.velocityZv.resize(b.vertexIndex.size()*b.nVertLevels);
	b.zTopVertex.resize(b.vertexIndex.size()*b.nVertLevels);

	
	for (int i=0;i<b.nCells;i++){



		for(int j=0;j<6;j++){
			int nvgid = b.verticesOnCell[i*6+j]; // neighbor vertex global id
			
				// get local vertex id from global id
			int nvlid = b.vertexIndex[nvgid];

			int c0 = b.cellIndex[b.cellsOnVertex[nvlid*3]], c1 = b.cellIndex[b.cellsOnVertex[nvlid*3+1]], c2 = b.cellIndex[b.cellsOnVertex[nvlid*3+2]];

				// if (mpas_c.indexToCellID[i]==4463){
				// 	dprint("4463 found in %d, nei nvgid vlid,  c0g c1g c2g, c0 c1 c2,  %d %d, %d %d %d, %d %d %d", mpas_c.gid, nvgid, nvlid, mpas_c.cellsOnVertex[nvlid*3], mpas_c.cellsOnVertex[nvlid*3+1], mpas_c.cellsOnVertex[nvlid*3+2], c0, c1, c2);
				// }



			double X[3][3] = {
				{b.xyzCell[c0*3], b.xyzCell[c0*3+1], b.xyzCell[c0*3+2]},
				{b.xyzCell[c1*3], b.xyzCell[c1*3+1], b.xyzCell[c1*3+2]},
				{b.xyzCell[c2*3], b.xyzCell[c2*3+1], b.xyzCell[c2*3+2]}
			};

			double P[3] = {b.xVertex[nvlid], b.yVertex[nvlid], b.zVertex[nvlid]};

			double lambda[3];
			barycentric_point2triangle(X[0], X[1], X[2], P, lambda);


			for (int curVertLevel=0;curVertLevel<b.nVertLevels; curVertLevel++){

				// mpas_c.zTopVertex[nvlid*mpas_c.nVertLevels + curVertLevel] = lambda[0] * mpas_c.zTop[c0*nVertLevels+curVertLevel] + lambda[1] * mpas_c.zTop[c1*nVertLevels+curVertLevel] + lambda[2] * mpas_c.zTop[c2*nVertLevels+curVertLevel];

				b.zTopVertex[nvlid*b.nVertLevels + curVertLevel] = lambda[0] * b.zTop[c0*nVertLevels+curVertLevel] + lambda[1] * b.zTop[c1*nVertLevels+curVertLevel] + lambda[2] * b.zTop[c2*nVertLevels+curVertLevel];


				if (b.cellsOnVertex[nvlid*3]==0 || b.cellsOnVertex[nvlid*3+1] ==0 || b.cellsOnVertex[nvlid*3+2] == 0){

					// mpas_c.velocityXv[nvlid*mpas_c.nVertLevels+curVertLevel] = 0;
					// mpas_c.velocityYv[nvlid*mpas_c.nVertLevels+curVertLevel] = 0;
					// mpas_c.velocityZv[nvlid*mpas_c.nVertLevels+curVertLevel] = 0;

					b.velocityXv[nvlid*b.nVertLevels+curVertLevel] = 0;
					b.velocityYv[nvlid*b.nVertLevels+curVertLevel] = 0;
					b.velocityZv[nvlid*b.nVertLevels+curVertLevel] = 0;
					
				}else{

					for (int curVertLevel=0;curVertLevel<b.nVertLevels; curVertLevel++){



						// mpas_c.velocityXv[nvlid*mpas_c.nVertLevels+curVertLevel] = lambda[0] * mpas_c.velocityX[c0*mpas_c.nVertLevels+curVertLevel] + lambda[1] * mpas_c.velocityX[c1*mpas_c.nVertLevels+curVertLevel] + lambda[2] * mpas_c.velocityX[c2*mpas_c.nVertLevels+curVertLevel];
						// mpas_c.velocityYv[nvlid*mpas_c.nVertLevels+curVertLevel] = lambda[0] * mpas_c.velocityY[c0*mpas_c.nVertLevels+curVertLevel] + lambda[1] * mpas_c.velocityY[c1*mpas_c.nVertLevels+curVertLevel] + lambda[2] * mpas_c.velocityY[c2*mpas_c.nVertLevels+curVertLevel];
						// mpas_c.velocityZv[nvlid*mpas_c.nVertLevels+curVertLevel] = lambda[0] * mpas_c.velocityZ[c0*mpas_c.nVertLevels+curVertLevel] + lambda[1] * mpas_c.velocityZ[c1*mpas_c.nVertLevels+curVertLevel] + lambda[2] * mpas_c.velocityZ[c2*mpas_c.nVertLevels+curVertLevel];
						
						b.velocityXv[nvlid*b.nVertLevels+curVertLevel] = lambda[0] * b.velocityX[c0*b.nVertLevels+curVertLevel] + lambda[1] * b.velocityX[c1*b.nVertLevels+curVertLevel] + lambda[2] * b.velocityX[c2*b.nVertLevels+curVertLevel];
						b.velocityYv[nvlid*b.nVertLevels+curVertLevel] = lambda[0] * b.velocityY[c0*b.nVertLevels+curVertLevel] + lambda[1] * b.velocityY[c1*b.nVertLevels+curVertLevel] + lambda[2] * b.velocityY[c2*b.nVertLevels+curVertLevel];
						b.velocityZv[nvlid*b.nVertLevels+curVertLevel] = lambda[0] * b.velocityZ[c0*b.nVertLevels+curVertLevel] + lambda[1] * b.velocityZ[c1*b.nVertLevels+curVertLevel] + lambda[2] * b.velocityZ[c2*b.nVertLevels+curVertLevel];	

					}


				}
				
			}
			

		}


	}



}

void populate_block( PBlock* b,
	int gid,
	std::vector<double> &data_bar, diy::mpi::communicator &world){

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




	}

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

	// mpas_c.indexToCellID = b->indexToCellID;
	// mpas_c.xCells = b->xCell;
	// mpas_c.yCells = b->yCell;
	// mpas_c.zCells = b->zCell;
	// mpas_c.velocityX = b->velocityX;
	// mpas_c.velocityY = b->velocityY;
	// mpas_c.velocityZ = b->velocityZ;
	// mpas_c.indexToVertexID = b->indexToVertexID;
	// mpas_c.xVertex = b->xVertex;
	// mpas_c.yVertex = b->yVertex;
	// mpas_c.zVertex = b->zVertex;
	// mpas_c.zTop = b->zTop;
	// mpas_c.cellsOnVertex = b->cellsOnVertex;
	// mpas_c.verticesOnCell = b->verticesOnCell;

	// mpas_c.nCells = mpas_c.xCells.size();
	// mpas_c.nVertices = mpas_c.xVertex.size();

	b->nCells = b->xCell.size();
	b->nVertices = b->xVertex.size();

	// mpas_c.gid = b->gid;




	// C_local.resize(3, nCells);
	// for (int i=0;i<nCells;i++){
	// 	C_local(0,i) = mpas_c.xCells[i];
	// 	C_local(1,i) = mpas_c.yCells[i];
	// 	C_local(2,i) = mpas_c.zCells[i];
	// }

	// mpas_c.nns_cells = Nabo::NNSearchD::createKDTreeLinearHeap(C_local);

	// M_local = new Eigen::MatrixXd(3, nVertices);
	// for (int i=0;i<nVertices;i++){
	// 	(*M_local)(0,i) = b->xVertex[i];
	// 	(*M_local)(1,i) = b->yVertex[i];
	// 	(*M_local)(2,i) = b->zVertex[i];
	// }


	// mpas_c.nns = Nabo::NNSearchD::createKDTreeLinearHeap(*M_local);

	
	// mpas_c.nVertLevels = 100; // FIX_ME: Remove the hard coding
	b->nVertLevels = 100; // FIX_ME: Remove the hard coding

	

	// mapping from global to local vertex id
	for (int i=0; i<b->nVertices; i++) {
		// mpas_c.vertexIndex[mpas_c.indexToVertexID[i]] = i;
		// mpas_c.block_vert_ids.insert(mpas_c.indexToVertexID[i]);

		b->vertexIndex[b->indexToVertexID[i]] = i;
		b->block_vert_ids.insert(b->indexToVertexID[i]);
		// fprintf(stderr, "%d, %d\n", i, indexToVertexID[i]);
	}


	// derive velocity on cell verticies 
	// mpas_c.velocityXv.resize(mpas_c.nVertices*mpas_c.nVertLevels);
	// mpas_c.velocityYv.resize(mpas_c.nVertices*mpas_c.nVertLevels);
	// mpas_c.velocityZv.resize(mpas_c.nVertices*mpas_c.nVertLevels);
	// mpas_c.zTopVertex.resize(mpas_c.nVertices*mpas_c.nVertLevels);
	// mpas_c.zTopVertexNorm.resize(mpas_c.nVertices*mpas_c.nVertLevels);
	// mpas_c.xyzCell.resize(mpas_c.nCells*3);

	b->velocityXv.resize(b->nVertices*b->nVertLevels);
	b->velocityYv.resize(b->nVertices*b->nVertLevels);
	b->velocityZv.resize(b->nVertices*b->nVertLevels);
	b->zTopVertex.resize(b->nVertices*b->nVertLevels);
	b->xyzCell.resize(b->nCells*3);


	for (int i=0; i<b->nCells; i++){
		// mpas_c.xyzCell[i*3] = mpas_c.xCells[i];
		b->xyzCell[i*3] = b->xCell[i];
	} 

	for (int i=0; i<b->nCells; i++){ 
		// mpas_c.xyzCell[i*3+1] = mpas_c.yCells[i];
		b->xyzCell[i*3+1] = b->yCell[i];

	}

	for (int i=0; i<b->nCells; i++){ 
		// mpas_c.xyzCell[i*3+2] = mpas_c.zCells[i];
		b->xyzCell[i*3+2] = b->zCell[i];
	}

	int max_cell_on_vert = 0;



}



void enq_halo_info( int cur_bid, int nprocs, 
	PBlock*                              b,
	const diy::Master::ProxyWithLink&   cp,
	const diy::Assigner& assigner){

	std::vector<std::set<int>> ghost_cell_ids; // 

	// std::vector<ghost_req> ghost_cell_reqs; // 
	
	

	ghost_cell_ids.resize(nprocs);
	
	// ghost_cell_reqs.resize(nprocs);




	// populate vector of required cell ids

	for (int i=0; i<b->nCells; i++){
		int cur_cell_1gid = b->indexToCellID[i];

		for (int j=0; j<b->cell_g_neighbors[cur_cell_1gid].size(); j++){
			
			// check if the current neighbor is in same partition/block
			int nbr_cell_1gid = b->cell_g_neighbors[cur_cell_1gid][j];

			int cell_bid = b->cgid_to_bid[nbr_cell_1gid];

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

void process_halo_req(PBlock* b, const diy::Master::ProxyWithLink& cp, const diy::Assigner& assigner){

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


				int cellIdxLocal = b->cellIndex[ghost_cell_ids[j]];

				// include this cell info in ghost point
				gp.cgid = ghost_cell_ids[j];
				int lgid = b->cellIndex[gp.cgid];
				gp.cx = b->xyzCell[3*lgid];
				gp.cy = b->xyzCell[3*lgid+1];
				gp.cz = b->xyzCell[3*lgid+2];



				// gp.vel_cx.insert(gp.vel_cx.end(), &mpas_c.velocityX[lgid*mpas_c.nVertLevels], &mpas_c.velocityX[(1+lgid)*mpas_c.nVertLevels]);
				// gp.vel_cy.insert(gp.vel_cy.end(), &mpas_c.velocityY[lgid*mpas_c.nVertLevels], &mpas_c.velocityY[(1+lgid)*mpas_c.nVertLevels]);
				// gp.vel_cz.insert(gp.vel_cz.end(), &mpas_c.velocityZ[lgid*mpas_c.nVertLevels], &mpas_c.velocityZ[(1+lgid)*mpas_c.nVertLevels]);
				for (int k=0;k<b->nVertLevels;k++){
					gp.vel_cx[k] = b->velocityX[lgid*b->nVertLevels + k];
					gp.vel_cy[k] = b->velocityY[lgid*b->nVertLevels + k];
					gp.vel_cz[k] = b->velocityZ[lgid*b->nVertLevels + k];
					gp.zTop[k] = b->zTop[lgid*b->nVertLevels + k];
				}



				int ngvid=-1; // neighbor global vertex id
				// get neighboring vertices
				for (int k=0;k<6;k++){
					// TODO: investigate crash on setting following to -1
					gp.vert_gids[k] = -1; // initially assume vertex not present
					
					
					ngvid = b->verticesOnCell[lgid*6+k];


					// if (b->gid==1){
					// 	fprintf(stderr, "spacking vertex %d %d\n", gp.cgid, ngvid);

					// 	} 
					

					// if vertex present in partition. include this vertex's info in ghost point.
					if (b->vertexIndex.find(ngvid)!=b->vertexIndex.end()){
						


						
						int nlvid = b->vertexIndex[ngvid]; // get local vid

						gp.vert_gids[k] = ngvid; // store global vertex id
						gp.vert_x[k] = b->xVertex[nlvid];
						gp.vert_y[k] = b->yVertex[nlvid];
						gp.vert_z[k] = b->zVertex[nlvid];

						// populate vertex neighbor cell-ids
						gp.cellsOnVertex[k*3] = b->cellsOnVertex[nlvid*3];
						gp.cellsOnVertex[k*3+1] = b->cellsOnVertex[nlvid*3+1];
						gp.cellsOnVertex[k*3+2] = b->cellsOnVertex[nlvid*3+2];

						// if (ngvid==14425){
						// 		fprintf(stderr, "found inend1 14425, %d %d %f\n", k, ngvid, gp.vert_x[k]);
						// 		}

					}

				}

				// push gp to ghost_info
				ghost_info.push_back(gp);


			}


			// queue cell info for exchange
			int dest_proc           = assigner.rank(dest_gid);
			diy::BlockID dest_block = {dest_gid, dest_proc};
			cp.enqueue(dest_block, ghost_info);


		}

		dprint("end of process process_halo_req in gid %d", b->gid);

	}

	void deq_halo_info(PBlock* b, const diy::Master::ProxyWithLink& cp){
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
					// mpas_c.xyzCell.push_back(ghost_info[j].cx);
					// mpas_c.xyzCell.push_back(ghost_info[j].cy);
					// mpas_c.xyzCell.push_back(ghost_info[j].cz);
					// std::copy(ghost_info[j].vel_cx, ghost_info[j].vel_cx+mpas_c.nVertLevels, std::back_inserter(mpas_c.velocityX));
					// std::copy(ghost_info[j].vel_cy, ghost_info[j].vel_cy+mpas_c.nVertLevels, std::back_inserter(mpas_c.velocityY));
					// std::copy(ghost_info[j].vel_cz, ghost_info[j].vel_cz+mpas_c.nVertLevels, std::back_inserter(mpas_c.velocityZ));
					// std::copy(ghost_info[j].zTop, ghost_info[j].zTop+mpas_c.nVertLevels, std::back_inserter(mpas_c.zTop));

					b->xyzCell.push_back(ghost_info[j].cx);
					b->xyzCell.push_back(ghost_info[j].cy);
					b->xyzCell.push_back(ghost_info[j].cz);
					std::copy(ghost_info[j].vel_cx, ghost_info[j].vel_cx+b->nVertLevels, std::back_inserter(b->velocityX));
					std::copy(ghost_info[j].vel_cy, ghost_info[j].vel_cy+b->nVertLevels, std::back_inserter(b->velocityY));
					std::copy(ghost_info[j].vel_cz, ghost_info[j].vel_cz+b->nVertLevels, std::back_inserter(b->velocityZ));
					std::copy(ghost_info[j].zTop, ghost_info[j].zTop+b->nVertLevels, std::back_inserter(b->zTop));

				// update indexToCellID, cellIndex
					// mpas_c.indexToCellID.push_back(ghost_info[j].cgid);
					// mpas_c.cellIndex[ghost_info[j].cgid] = mpas_c.indexToCellID.size()-1;

					b->indexToCellID.push_back(ghost_info[j].cgid);
					b->cellIndex[ghost_info[j].cgid] = b->indexToCellID.size()-1;


				// update ghost cell ids into std::set
					for (int k=0;k<6;k++){
						if (ghost_info[j].vert_gids[k]>-1){
						// process associated vertex

						// vertexIndex
							int vgid = ghost_info[j].vert_gids[k]; 

							if (b->vertexIndex.count(vgid) == 0){

								
								
								// mpas_c.vertexIndex[vgid] = (int) mpas_c.xVertex.size();
								b->vertexIndex[vgid] = (int) b->xVertex.size();

								// if (ghost_info[j].vert_gids[k]==14425){
								// 	dprint("found inend 14425, %d %d %f %f %f %d %d %d", mpas_c.vertexIndex[vgid], k, ghost_info[j].vert_x[k], ghost_info[j].vert_y[k], ghost_info[j].vert_z[k], ghost_info[j].cellsOnVertex[k*3], ghost_info[j].cellsOnVertex[k*3+1], ghost_info[j].cellsOnVertex[k*3]+2);
								// }

								// // vertex positions
								// mpas_c.xVertex.push_back(ghost_info[j].vert_x[k]);
								// mpas_c.yVertex.push_back(ghost_info[j].vert_y[k]);
								// mpas_c.zVertex.push_back(ghost_info[j].vert_z[k]);

								b->xVertex.push_back(ghost_info[j].vert_x[k]);
								b->yVertex.push_back(ghost_info[j].vert_y[k]);
								b->zVertex.push_back(ghost_info[j].vert_z[k]);

								// // cellOnVertex?*not needed for pure ghost vertices. needed for shared ghost. keep for all.
								// mpas_c.cellsOnVertex.push_back(ghost_info[j].cellsOnVertex[k*3]);
								// mpas_c.cellsOnVertex.push_back(ghost_info[j].cellsOnVertex[k*3+1]);
								// mpas_c.cellsOnVertex.push_back(ghost_info[j].cellsOnVertex[k*3+2]);

								b->cellsOnVertex.push_back(ghost_info[j].cellsOnVertex[k*3]);
								b->cellsOnVertex.push_back(ghost_info[j].cellsOnVertex[k*3+1]);
								b->cellsOnVertex.push_back(ghost_info[j].cellsOnVertex[k*3+2]);

							} 


						}
					}

				}


			}
		}

	// compute velocities at shared vertices after halo exchange. Only native and shared.
	}

// consumer
	void con(Decaf* decaf)
	{      	
		int max_rounds = 3;
	// int max_steps = 2000;

		int max_steps = 4000;
		std::string ip_file = "graph.topology";

		// mpaso mpas_c;
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
		

		// keep pointer to current block. use it for populating / updating only current block later
		b = new PBlock;
		b->generate_domain_decomposition_graph(ip_file, world.size());	
		b->read_cell_g_neighbors();
		pgid = gid;
		global_gid = gid;
		b->gid = gid; // TODO: remove other gid vars e.g. global_gid
		master.add(gid, b, link);    // add the current local block to the master
	}

	
	// mpas_g.generate_domain_decomposition_graph(ip_file, world.size());	
	// mpas_g.read_cell_g_neighbors();



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
				
				int frame_no = (int) data_bar[14];

				if (frame_no==1){

					// streamlines slines;
					b->compute_cellIndex(b->gid,world.size());

					// TODO  premature stopping if work completed before max_rounds 
					populate_block(b, pgid, data_bar, world);

					// identify and enqueue halo cells
					master.foreach([&](PBlock* b, const diy::Master::ProxyWithLink& cp)
						{ enq_halo_info(b->gid, world.size(), 
							b, cp, assigner); });
					
					//rexchange
					bool remote = true;
					master.exchange(remote);
					// process the halo request
					master.foreach([&](PBlock* b, const diy::Master::ProxyWithLink& cp){
						process_halo_req(b, cp, assigner);
					});
					// rexchange
					master.exchange(remote);

					//deq_halo_info
					master.foreach([&](PBlock* b, const diy::Master::ProxyWithLink& cp){
						deq_halo_info(b, cp);
					});

					compute_ghost_vertex_values(*b);

					int stop = (max_rounds ? max_rounds : 1);
					int incr = (max_rounds ? 1 : 0);
					for (int round = 0; round < stop; round += incr)
					{

						master.foreach([&](PBlock* b, const diy::Master::ProxyWithLink& cp)
						{

							parallel_streamlines( round, *b, max_steps, cp,  assigner, world);
						});
						bool remote = true;
						master.exchange(remote);	
					}

					//bio.write_cell_centers(pgid, world, b->xCell); //only testing
					bio.write_particle_traces(pgid, world, *b, max_steps);
				}	

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
bool get_curpos_vel_sl( PBlock &b, double *X, Eigen::Vector3d &c_vel, bool &in_global_domain){

	int K = 6;
	double radius = b.radius;
	double depth = std::sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]) - radius;



	// get nearest vertex neighbor ids using mpas_c
	Eigen::VectorXi nearest_idx(K);
	Eigen::VectorXd dists2(K);

	Eigen::VectorXd q(3);
	q<<X[0], X[1], X[2];

	

	// get nearest cell neighbor ids using mpas_c
	Eigen::VectorXi nearest_cell_idx(1);
	Eigen::VectorXd dists2_cell(1);
	b.nns_cells->knn(q, nearest_cell_idx,dists2_cell, 1,0, Nabo::NNSearchF::SORT_RESULTS);

	int clid = nearest_cell_idx[0]; // cell local id
	for (int i=0;i<K;i++){
		int gvid = b.verticesOnCell[clid*6+i];
		nearest_idx[i] = b.vertexIndex[gvid];
	}
	// dprint("clid %d %d %f %f %f", clid, mpas_c.indexToCellID[clid], q[0], q[1], q[2]);
	if (clid>b.nCells-1){
		return false; // out of partition
	}


	// interpolate vertically and get positions and velocity (in values)
	std::vector<double> values;
	values.resize(6*nearest_idx.size());


	double top_boun = b.zTop[b.nVertLevels*nearest_cell_idx[0]+0];

	double botom_boun = b.zTop[b.nVertLevels*nearest_cell_idx[0]+99];


	if (depth> top_boun || depth< botom_boun){
		// terminate flow line
		in_global_domain = false;
		dprint("vertical exit printing values: %f %f %f", top_boun, botom_boun, depth);
		return false; 
		
	}else{ 
		interpolate_vertically(b.nVertLevels, b.zTopVertex, nearest_idx, values, depth, b.xVertex, b.yVertex, b.zVertex, b.velocityXv, b.velocityYv, b.velocityZv, b.gid);

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
	int cgid = b.indexToCellID[clid];
	if (b.gid == 0 && cgid == 6){
		dprint("nearest_cell_idx for 0 clid cgid: %d %d, %d %d %d %d %d %d, %f %f %f", clid, cgid, nearest_idx[0], nearest_idx[1], nearest_idx[2], nearest_idx[3], nearest_idx[4], nearest_idx[5], c_vel[0], c_vel[1], c_vel[2] );
	}


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


bool rk1_mpas(
	double *X,
	double h,
	double *Y,
	bool &in_global_domain, 
	PBlock &b)
{

	
	Eigen::Vector3d c_vel;

	if (get_curpos_vel_sl( b, X, c_vel, in_global_domain)){	

		// advect
		Y[0] = X[0] + h*c_vel[0];	
		Y[1] = X[1] + h*c_vel[1];
		Y[2] = X[2] + h*c_vel[2];

		Eigen::MatrixXd R = project_on_sphere(X, Y, Y);


		return true;

	}

	return false;
	

}

void parallel_streamlines( int round, PBlock &b, int max_steps, const diy::Master::ProxyWithLink&   cp, const diy::Assigner& assigner, diy::mpi::communicator &world){

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
		int n_seed_verts = b.nVertices/skipval;
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


				dc_c = b.radius + b.zTopVertex[b.nVertLevels*vert_id+seed_level_ids[j]];
				dc_t = b.radius + b.zTopVertex[b.nVertLevels*vert_id+0];
				dc_b = b.radius + b.zTopVertex[b.nVertLevels*vert_id+99];




				dc_base = std::sqrt(b.xVertex[vert_id]*b.xVertex[vert_id] + b.yVertex[vert_id]*b.yVertex[vert_id] + b.zVertex[vert_id]*b.zVertex[vert_id]);
				ratio = dc_c/dc_base;

				EndPt p;
				p.pid = b.init;
				p.sid = b.init;
				p.step = 0;
				p.gpid = b.global_start_ids[b.gid] + p.pid;
				p[0] = b.xVertex[vert_id] * ratio;
				p[1] = b.yVertex[vert_id] * ratio;
				p[2] = b.zVertex[vert_id] * ratio;

				dprint("vlid gvid bid %d %d %d", vert_id, b.indexToVertexID[vert_id], b.gid);

				particles.push_back(p);
				b.init++;



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
		
		
		while(rk1_mpas(cur_p.coords, 100000, next_p.coords, in_global_domain, b)){	
		// while(rk1_mpas(mpas_g, mpas_c, cur_p.coords, 10000, next_p.coords, in_global_domain)){	
			s.pts.push_back(next_p);
			cur_p = next_p;


			// TODO: is the nan check really needed (deal with nan earlier skip_val=3000)?
			if (s.start_step + s.pts.size() >= max_steps || std::isnan(next_p.coords[0]) || std::isnan(next_p.coords[1]) || std::isnan(next_p.coords[2])) 
			{
				dprint("Reached max steps");
				finished = true;
				break;

			} 


		}

		b.global_trace_sizes[s.gpid] = s.start_step + s.pts.size();


		if (in_global_domain==false)
			finished = true;


		b.segments.push_back(s);	

		if (finished){
			b.done++;

			dprint("finished particle in %d", b.gid);
		}
		else{ 			// package unfinished endpoint for sending 
			dprint("enqueuing particle");
			EndPt out_pt(s); // TO_DO: Also update the start_time here
			// int dest_gid = mpas_c.get_bid_for_pt(next_p.coords);		
			int dest_gid = b.get_bid_for_pt(next_p.coords);		
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



