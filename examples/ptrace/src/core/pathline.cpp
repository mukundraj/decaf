#include "pathline.h"
#include "mpas_io.h"
#include "misc.h"
#include <Eigen/Dense>
#include "geometry_utils.h"
#include "block.h"
#include "advect.h"
#include <unordered_set>

using namespace Eigen;

bool pathline::in_global_domain(const Pt& p, double cx, double cy, double cz){

	// if (pow(-6352390-p.coords[0],2)+pow(-447879-p.coords[1],2)+(pow(197714-p.coords[2],2))>pow(5524190.98744,2))
	// 	return false;
	// else
	// 	return true;

	// if (pow(cx-p.coords[0],2)+pow(cy-p.coords[1],2)+(cz-p.coords[2],2)>pow(65000,2) || pow(p.coords[0],2)+pow(p.coords[1],2)+(p.coords[2],2) > pow(radius + 10, 2))
	// if (pow(cx-p.coords[0],2)+pow(cy-p.coords[1],2)+(cz-p.coords[2],2)>pow(65000,2))
	if (pow(cx-p.coords[0],2)+pow(cy-p.coords[1],2)+(cz-p.coords[2],2)>pow(130000,2))
		return false;
	else
		return true;
}

bool pathline::in_local_domain (const block *b, const Pt& p, int &iCell, int round){

	int nCellVertices;
	Array3d xSubStep;
	xSubStep(0) = p.coords[0];
	xSubStep(1) = p.coords[1];
	xSubStep(2) = p.coords[2];
	get_validated_cell_id(*b, xSubStep, iCell, nCellVertices);
	// dprint("bgid %d iCell %d b->in_partition[iCell] %d round %d  b->gcIdxToGid[iCell] %d", b->gid, iCell, b->in_partition[iCell], round,  b->gcIdxToGid[iCell]);
	
	// if (b->gcIdxToGid[iCell] == b->gid)
	if (b->in_partition[iCell] == round)
		return true;
	else
		return false;
}

bool pathline::compute_streamlines(block *b, mpas_io &mpas1, const diy::Master::ProxyWithLink &cp, const diy::Assigner &assigner, int prediction, size_t &nsteps, std::vector<EndPt> &particles_hold, int skip_rate){

	radius = b->radius;

	Array2d kWeightK = {0, 0.5}, kWeightT = {0, 0.5};
	Array2d kWeightKVert = kWeightK, kWeightTVert = kWeightT;

	ArrayXXd kWeightX(3, 2);
	kWeightX << 0, 1,
		0, 1,
		0, 1;
	int subStepOrder = 2;
	Array2d kWeightXVert = kWeightX.block<1, 2>(0, 0); // see https://eigen.tuxfamily.org/dox/group__TutorialBlockOperations.html
	// Array2d kWeightXVert{0,1};

	// local variables
	Array3d particleVelocity, xSubStep, diffSubStep, diffParticlePosition;
	double tSubStep, diffParticlePositionVert;
	double zSubStep, zLevelParticle, diffSubStepVert;
	int iCell = -1, nCellVertices, iLevel;
	double particleVelocityVert;

	size_t nCellsnVertLeves = b->nCells * b->nVertLevels;
	dprint("b->particles.size() %ld, nCells %ld", b->particles.size(), b->nCells);
	int iCell_prev;

	std::vector<EndPt> particles_finished;

		for (int pid = 0; pid < b->particles.size(); pid++)
		{
			// create segment
			Segment seg;
			seg.pid = b->particles[pid].pid;
			seg.nsteps = b->particles[pid].nsteps;

			Array3d particlePosition = {
				b->particles[pid][0],
				b->particles[pid][1],
				b->particles[pid][2]};
			int &cur_nsteps = b->particles[pid].nsteps;
			// read zLevelParticle
			zLevelParticle = b->particles[pid].zLevelParticle; // particles_zLevel[0];
			// iCell = -1;											   //mpas1->particles[pid].glCellIdx;

			iCell = b->particles[pid].glCellIdx;

			bool finished = false;

			if (b->particles[pid].pid==1000)
							dprint("next particle gid %d pid %d, cur_nsteps %d, nSteps %f | iCell %d, (%f %f %f),zLevP %f", b->gid, b->particles[pid].pid, cur_nsteps, nSteps, iCell, particlePosition(0), particlePosition(1), particlePosition(2), zLevelParticle);

			
			// if (b->particles[pid].pid!=5470)
			// 	{
			// 		finished = true;
			// 		continue;
			// 	}
			// dprint("starting %d | %f %f %f | %f", b->particles[pid].pid, b->particles[pid][0], b->particles[pid][1], b->particles[pid][2], zLevelParticle);

			for (int timeStep = 0; timeStep < nSteps; timeStep++)
			{
					ArrayXXd kCoeff = ArrayXXd::Zero(3, subStepOrder + 1);
					ArrayXXd kCoeffVert = ArrayXXd::Zero(1, subStepOrder + 1);
					for (int subStep = 0; subStep < subStepOrder; subStep++)
					{	
						xSubStep = particlePosition;
						diffSubStep = kWeightK(subStep) * kCoeff.col(subStep);
						if (kWeightK(subStep) != 0)
						{
							particle_horizontal_movement(*b, xSubStep, diffSubStep);
						}
						
						// vertical
						zSubStep = zLevelParticle;
						diffSubStepVert = kWeightKVert(subStep) * kCoeffVert(subStep);


						if (kWeightKVert(subStep) != 0)
						{
							zSubStep = zSubStep + diffSubStepVert;
						}
						// dprint("before get_validated..");
						// get new time step (tm = (timestep-1)*dt)
						tSubStep = (timeStep + kWeightT(subStep)) * dt;
						get_validated_cell_id(*b, xSubStep, iCell, nCellVertices);
						// iLevel = get_vertical_id(b->maxLevelCell[iCell], zSubStep, &zMid_cur[iCell * b->nVertLevels]);
						iLevel = get_vertical_id(b->maxLevelCell[iCell], zSubStep, &b->zMid[iCell * b->nVertLevels]);
						if (subStep==0)
							iCell_prev = iCell;

						// dprint("iCell   %d iLevel %d, nSteps %f", iCell, iLevel, nSteps);
						// dprint("tStep %d, subStep %d, pid %d, iCell %d, iLevel %d, maxLevelCell %d, zSubStep %f, zMid[iLevel] %f, zMid[iLevel+1] %f", timeStep, subStep, b->particles[pid].pid, iCell, iLevel, b->maxLevelCell[iCell], zSubStep, zMid_cur[iCell * b->nVertLevels+iLevel], zMid_cur[iCell * b->nVertLevels+iLevel+1]);

						int timeInterpOrder = 2;
						double timeCoeff[2];
						timeCoeff[0] = tSubStep / dtSim;
						timeCoeff[1] = 1.0 - timeCoeff[0];
						velocity_time_interpolation(timeInterpOrder,
											timeCoeff,
											iCell,
											iLevel,
											nCellVertices,
											zSubStep,
											b,
											mpas1,
											xSubStep,
											particleVelocity,
											particleVelocityVert);
						// dprint("vel %d %d| %f %f %f", iCell, iLevel, particleVelocity(0), particleVelocity(1), particleVelocity(2));
						//!!!!!!!!!! FORM INTEGRATION WEIGHTS kj !!!!!!!!!!
						kCoeff.col(subStep + 1) = dt * particleVelocity;
						kCoeffVert(subStep + 1) = dt * particleVelocityVert;

					} // subStep loop ends

						// update particle positions
						diffParticlePosition << 0, 0, 0;
						diffParticlePositionVert = 0;

						for (int subStep = 0; subStep < subStepOrder; subStep++)
						{

							//first complete particle integration
							diffParticlePosition = diffParticlePosition + kWeightX.col(subStep) * kCoeff.col(subStep + 1);
							diffParticlePositionVert = diffParticlePositionVert + kWeightXVert(subStep) * kCoeffVert(subStep + 1);
						}



						//now, make sure particle position is still on same spherical shell as before
						particle_horizontal_movement(*b, particlePosition, diffParticlePosition);

						//now can do any vertical movements independent of the horizontal movement
						zLevelParticle = zLevelParticle + diffParticlePositionVert;

						// dprint(" timestep03 %d, particlePosition (%.7f %.7f %.7f), zLevelParticle %.7f", timeStep, particlePosition[0], particlePosition[1],particlePosition[2], zLevelParticle);

						// dprint("pp (%d %d) %f %f %f | %f", iCell, iLevel, particlePosition[0], particlePosition[1], particlePosition[2], zLevelParticle);

						// update segment
						Pt p;
						p.coords[0] = particlePosition[0]; 
						p.coords[1] = particlePosition[1]; 
						p.coords[2] = particlePosition[2];
						// p.zLevels.push_back(zLevelParticle);
						// if (cur_p_step %100 == 0)
							seg.pts.push_back(p);

						cur_nsteps ++;
						nsteps++;
						


						// if predicting, add copy coordinates to EndPt and add to b->particles_hold
						if (prediction == true && cur_nsteps % skip_rate == 0)
						{   
						    EndPt way_pt;
						    way_pt[0] = b->xCell[iCell];
						    way_pt[1] = b->yCell[iCell];
						    way_pt[2] = b->zCell[iCell];
						    way_pt.predonly = 1;
						    particles_hold.push_back(way_pt);
						}

						if (iCell_prev == 52571 || iCell_prev == 135779 || iCell_prev == 151446 || iCell_prev == 172732 || iCell_prev == 189476 || iCell_prev == 115789 || iCell_prev == 201643 || iCell_prev == 158103 || iCell_prev == 191753 || iCell_prev == 30044 || iCell_prev == 51769 || iCell_prev == 126417 || iCell_prev == 54758 || iCell_prev == 76847)
						{
							// dprint("particle stopped gid %d, pid %d steps %d", b->gid, b->particles[pid].pid, b->particles[pid].nsteps);
							finished = true;
							break;
						}

						// if (b->particles[pid].pid==1000)
						// 	dprint("next particle gid %d pid %d, cur_nsteps %d, nSteps %f | iCell %d, (%f %f %f), zLevP %f", b->gid, b->particles[pid].pid, cur_nsteps, nSteps, iCell, particlePosition(0), particlePosition(1), particlePosition(2), zLevelParticle);

						

						// check if particle inside global domain and handle
						if (!in_global_domain(p, b->xCell[iCell], b->yCell[iCell], b->zCell[iCell] ) || cur_nsteps >= nSteps)
						// if (cur_nsteps >= nSteps || iCell == 31823)
						{
							finished = true;
						
							
						}

						if (finished==true){
							// if (prediction==true)

							// mpas1.debugarray[b->particles[pid].pid/1000] += cur_nsteps;
							// dprint("finned %d in %d, pred %d | %d", b->particles[pid].pid, cur_nsteps, prediction, mpas1.debugarray[0]);
							break;
						}
						
						// dprint("nsteps %d, finished %d", nsteps, finished);
						// check if inside local domain and handle
						int round = 1;
						// dprint("iCell %d", iCell);
						if (!in_local_domain(b , p, iCell, round)){
				            if (b->particles[pid].pid==1000)
							dprint("jumped!! in %d, cur_nsteps %d, pid %d", b->gid, cur_nsteps, b->particles[pid].pid);
							break;
						}


			} // timeStep loop ends


			// if unfinished, enqueue 
			if (finished == false){
				int dest_gid = b->gcIdToGid[iCell+1];

					int dest_proc = assigner.rank(dest_gid);
					diy::BlockID dest_block = {dest_gid, dest_proc};
					// dprint("sou %d des %d, iCell (%d %d), cur_nsteps %d, pid %d, b->gcIdxToGid[iCell] %d, b->in_partition[iCell] %d",  b->gid, dest_gid, iCell_prev, iCell, cur_nsteps, b->particles[pid].pid, b->gcIdxToGid[iCell], b->in_partition[iCell]);
					EndPt pt;
					pt.pid = b->particles[pid].pid; // Needs modification of diy code to be effective
					pt[0] = particlePosition[0];  pt[1] = particlePosition[1];  pt[2] = particlePosition[2];
					pt.nsteps = b->particles[pid].nsteps;
					pt.zLevelParticle = zLevelParticle;
					// mpas1->particles[pid].glCellIdx = iCell; 
					pt.glCellIdx = iCell;

					diy::Link *l = static_cast<diy::Link *>(cp.link());
					bool jumped = true;
					for (size_t i = 0; i < l->size(); ++i)
					{
						int nbr_gid = l->target(i).gid;
						if (nbr_gid == dest_gid){
							jumped = false;
							cp.enqueue(dest_block, pt);
							break;
						}

					}
					if (jumped)
					{
						dprint("alert! jump gid %d dest_gid %d, pid %d, prev iCell %d, iCell %d", b->gid, dest_gid, b->particles[pid].pid, b->particles[pid].glCellIdx, iCell);
					}
			}

			// push segment to block segment vector
			// dprint("segsize %ld", seg.pts.size());
			// dprint ("seg.nsteps %d", seg.nsteps);
			b->segments.push_back(seg);
			// dprint("numsegs %ld", b->segments.size());

		} // particle loop

		// dprint("b->segments size %ld , gid %d", b->segments.size(), cp.gid());
	

	if (b->particles_store.size() > 0 ){
		// dprint("moving from store in %d", b->gid);
		b->particles = std::move(b->particles_store);
		// dprint("finishedt a callback in gid %d, %ld %ld", b->gid,  b->particles.size(), b->particles_store.size());
		return false;
	}
	else {
		// b->particle_store.clear();
		// dprint("finishedf a callback in gid %d, %ld %ld", b->gid,  b->particles.size(), b->particles_store.size());
		b->particles.clear(); // clearing particles if no particles came in
		return true;
	}
		
	

}


void pathline::cell_to_vertex_interpolation(block *b, std::vector<int> &local_gcIds){

	std::unordered_set<int> completed_verts;

	// dprint("verticesOnCell size %ld", b->verticesOnCell.size());

	

	for (size_t i=0; i<local_gcIds.size(); i++){
		int cellIdx = local_gcIds[i]-1;

		

		// loop over vertices on cellId
		for (int j=0; j<b->maxEdges; j++){
			int vid = b->verticesOnCell[b->maxEdges*cellIdx + j];

			// if (cellIdx == 5470)
			// 		dprint("FOUND 5470, vid %d | %d | size %ld", vid, completed_verts.find(vid)==completed_verts.end(), completed_verts.size());
			
			// if vertId in completed_verts, skip remaining in this iteration
			if (completed_verts.find(vid)==completed_verts.end()){
				completed_verts.insert(vid);
				int aVertex = vid - 1;
				int flag=0;
				// if (cellIdx == 5470){
				// 	dprint("FOUND 5470, aVertex %d", aVertex);
				// 	flag = 1;
				// }
					
				
				const int vertexDegree = 3;
				
				// get pointVertex: pos of cells surrounding vertex
				double pointVertex[vertexDegree][3];
				for (int aCell=0; aCell < vertexDegree; aCell++){

					int cellIdx = b->cellsOnVertex[vertexDegree*aVertex + aCell] - 1;
					// pointVertex[0][aCell] = b->xCell[b->cellsOnVertex[cellIdx ]];
					// pointVertex[1][aCell] = b->yCell[b->cellsOnVertex[cellIdx ]];
					// pointVertex[2][aCell] = b->zCell[b->cellsOnVertex[cellIdx ]];
					pointVertex[aCell][0] = b->xCell[cellIdx ];
					pointVertex[aCell][1] = b->yCell[cellIdx ];
					pointVertex[aCell][2] = b->zCell[cellIdx ];

					// if (flag==1)
					// {
					// 	dprint("cellIdx %d, (%f %f %f) %f", cellIdx, b->xCell[cellIdx ], b->yCell[cellIdx ], b->zCell[cellIdx ], pointVertex[0][0]);
					// }

				}

				// get pointInterp: pos of vertex
				Array3d pointInterp;
				pointInterp(0) = b->xVertex[aVertex];
				pointInterp(1) = b->yVertex[aVertex];
				pointInterp(2) = b->zVertex[aVertex];

				// get interpolation constants (lambda)
				double lambda[vertexDegree];
				double areaB[vertexDegree];
				wachspress_coordinates(*b, vertexDegree, pointVertex, pointInterp, areaB, &lambda[0], flag);

				// if (b->gid == 0)
				// if (flag == 1)
				// dprint("lambda %f %f %f, pointInterp %f %f %f | (%f %f %f) (%f %f %f) (%f %f %f) ", lambda[0], lambda[1], lambda[2], pointInterp(0), pointInterp(1), pointInterp(2), pointVertex[0][0], pointVertex[1][0], pointVertex[2][0], pointVertex[0][1], pointVertex[1][1], pointVertex[2][1], pointVertex[0][2], pointVertex[1][2], pointVertex[2][2]);
				
				// interpolate and store
				for (int aLevel=0; aLevel<b->nVertLevels; aLevel++){
					if (b->boundaryVertex[b->nVertLevels*aVertex + aLevel] < 1){
						
						// get cell center velocities
						double ucCell[vertexDegree][3];

						for(int aCell=0; aCell < vertexDegree; aCell++){
							int cellIdx = b->cellsOnVertex[vertexDegree*aVertex + aCell] - 1;
							ucCell[aCell][0] = b->velocityX[cellIdx * b->nVertLevels + aLevel];
							ucCell[aCell][1] = b->velocityY[cellIdx * b->nVertLevels + aLevel];
							ucCell[aCell][2] = b->velocityZ[cellIdx * b->nVertLevels + aLevel];
							// if (cellIdx == 1000){
							// 	dprint("ucCell %f %f %f", ucCell[aCell][0], ucCell[aCell][1], ucCell[aCell][2]);
							// }
						}

						b->uVertexVelocity[b->nVertLevels*aVertex + aLevel] = wachspress_interpolate(lambda, ucCell, 0, vertexDegree);
						b->vVertexVelocity[b->nVertLevels*aVertex + aLevel] = wachspress_interpolate(lambda, ucCell, 1, vertexDegree);
						b->wVertexVelocity[b->nVertLevels*aVertex + aLevel] = wachspress_interpolate(lambda, ucCell, 2, vertexDegree);
					}
				}


			}

		}

			


	}

	
}

pathline::pathline(mpas_io &mpas1, double dtSim_in, double dtParticle_in)
{

	uVertexVelocities.resize(timeInterpOrder);
	vVertexVelocities.resize(timeInterpOrder);
	wVertexVelocities.resize(timeInterpOrder);

	dtSim = dtSim_in;
	dtParticle = dtParticle_in;

	// adjust time step for consistency with integer number of steps
	nSteps = ceil(dtSim / dtParticle);
	dt = dtSim / nSteps;
}

void pathline::compute_epoch(block *mpas1, int framenum)
{

// 	Array2d kWeightK = {0, 0.5}, kWeightT = {0, 0.5};
// 	Array2d kWeightKVert = kWeightK, kWeightTVert = kWeightT;

// 	ArrayXXd kWeightX(3, 2);
// 	kWeightX << 0, 1,
// 		0, 1,
// 		0, 1;
// 	int subStepOrder = 2;
// 	Array2d kWeightXVert = kWeightX.block<1, 2>(0, 0); // see https://eigen.tuxfamily.org/dox/group__TutorialBlockOperations.html
// 	// Array2d kWeightXVert{0,1};

// 	// local variables
// 	Array3d particleVelocity, xSubStep, diffSubStep, diffParticlePosition;
// 	double tSubStep, diffParticlePositionVert;
// 	double zSubStep, zLevelParticle, diffSubStepVert;
// 	int iCell = -1, nCellVertices, iLevel;
// 	double particleVelocityVert;

// 	size_t nCellsnVertLeves = mpas1->nCells * mpas1->nVertLevels;

// 	// dprint("mpas1->particles.size() %ld, nCells %ld", mpas1->particles.size(), mpas1->nCells);

// 	std::vector<EndPt> particles_finished;

// 	for (int pid = 0; pid < mpas1->particles.size(); pid++)
// 	{

// 		// dprint("pid %d", pid);
// 		// create segment
// 		Segment seg;
// 		seg.pid = mpas1->particles[pid].pid;

// 		// if (pid != 8 ) continue;
// 		// if (mpas1->gid != 1 ) continue;

// 		// Array3d particlePosition = {
// 		// 	particles_current[3*pid+0],
// 		// 	particles_current[3*pid+1],
// 		// 	particles_current[3*pid+2]};

// 		Array3d particlePosition = {
// 			mpas1->particles[pid][0],
// 			mpas1->particles[pid][1],
// 			mpas1->particles[pid][2]};

// 		// read zLevelParticle
// 		// zLevelParticle = mpas1->zLevelParticle[simStep*nParticles+pid];
// 		// zLevelParticle = zparticle_current[pid];
// 		zLevelParticle = mpas1->particles[pid].zLevelParticle; // particles_zLevel[0];
// 		// iCell = -1;											   //mpas1->particles[pid].glCellIdx;

		

// 		iCell = mpas1->particles[pid].glCellIdx;

// 		// dprint("pid %d zLevelParticle %f, particlePosition (%f %f %f)", pid, zLevelParticle, particlePosition[0], particlePosition[1], particlePosition[2]);

// 		// dprint("zLevelParticle %f", zLevelParticle);
// 		bool finished = false;
// 		for (int timeStep = 0; timeStep < nSteps; timeStep++)
// 		{

// 			ArrayXXd kCoeff = ArrayXXd::Zero(3, subStepOrder + 1);
// 			ArrayXXd kCoeffVert = ArrayXXd::Zero(1, subStepOrder + 1);

			
// 			for (int subStep = 0; subStep < subStepOrder; subStep++)
// 			{

// 				xSubStep = particlePosition;
// 				diffSubStep = kWeightK(subStep) * kCoeff.col(subStep);

// 				if (kWeightK(subStep) != 0)
// 				{
// 					particle_horizontal_movement(*mpas1, xSubStep, diffSubStep);
// 				}

// 				// vertical
// 				zSubStep = zLevelParticle;
// 				diffSubStepVert = kWeightKVert(subStep) * kCoeffVert(subStep);

// 				// dprint("timestep00 %d, subStep %d, pid %d, (%.7f %.7f %.7f), diff % f (%f %f %f) zLevel %.7f", timeStep, subStep, mpas1->particles[pid].pid, xSubStep(0), xSubStep(1), xSubStep(2), diffSubStepVert, diffSubStep(0), diffSubStep(1), diffSubStep(2), zLevelParticle);


// 				if (kWeightKVert(subStep) != 0)
// 				{
// 					zSubStep = zSubStep + diffSubStepVert;
// 				}

// 				// get new time step (tm = (timestep-1)*dt)
// 				tSubStep = (timeStep + kWeightT(subStep)) * dt;

// 				get_validated_cell_id(*mpas1, xSubStep, iCell, nCellVertices);

// 				// dprint("pid %d, Cell %d, gid %d", pid, iCell, mpas1->gid);
				

// 				// get cellIdx

// 				// iCell = mpas1->cellIndex[iCell];
// 				iCell = iCell;
				

// 				// dprint("iCell %d %f", iCell, mpas1->xCell[iCell]);

// 				// dprint("zMid_cur %ld, zTop_cur %ld", zMid_cur_size, zTop_cur_size);
// 				iLevel = get_vertical_id(mpas1->maxLevelCell[iCell], zSubStep, &zMid_cur[iCell * mpas1->nVertLevels]);
// 				// dprint("iLevel %d", iLevel);

// 				int timeInterpOrder = 2;
// 				double timeCoeff[2];
// 				timeCoeff[0] = tSubStep / dtSim;
// 				timeCoeff[1] = 1.0 - timeCoeff[0];
// 				// dprint("mid");
// 				// dprint("tSubStep %f, timeCoeff (%f %f)", tSubStep, timeCoeff[0], timeCoeff[1]);

// 				// velocity_time_interpolation
// 				// double v[3];
// 				// double p[3] = {xSubStep(0), xSubStep(1), xSubStep(2)};

// 				velocity_time_interpolation(timeInterpOrder,
// 											timeCoeff,
// 											iCell,
// 											iLevel,
// 											nCellVertices,
// 											zSubStep,
// 											b,
// 											*mpas1,
// 											xSubStep,
// 											particleVelocity,
// 											particleVelocityVert);

// 				// dprint(" %f (%f %f %f)", zLevelParticle, particleVelocity[0], particleVelocity[1],particleVelocity[2]); 
// 				// dprint("timestep01 %d,subStep %d, pid %d,(%f %f %f), velVert %.07e, particleVelocity (%.07e %.7ef %.7e), zLevel %f", timeStep, subStep, mpas1->particles[pid].pid, xSubStep(0), xSubStep(1), xSubStep(2), particleVelocityVert, particleVelocity(0), particleVelocity(1), particleVelocity(2), zLevelParticle);


// 				//!!!!!!!!!! FORM INTEGRATION WEIGHTS kj !!!!!!!!!!
// 				kCoeff.col(subStep + 1) = dt * particleVelocity;
// 				kCoeffVert(subStep + 1) = dt * particleVelocityVert;

// 				// dprint(" %f (%f )", zLevelParticle, kCoeff.col(subStep + 1)(0)); 
// 				iCell = mpas1->indexToCellID[iCell];
// 			} // substep loop ends

			

// 			// update particle positions
// 			diffParticlePosition << 0,0,0;
// 			diffParticlePositionVert = 0;

// 			for (int subStep = 0; subStep < subStepOrder; subStep++){

//               		//first complete particle integration
// 					diffParticlePosition = diffParticlePosition + kWeightX.col(subStep) * kCoeff.col(subStep+1);
// 					diffParticlePositionVert = diffParticlePositionVert + kWeightXVert(subStep) * kCoeffVert(subStep+1);
// 			}

// 			// dprint("timestep02 diffParticlePosition %.7e %.7e %.7e", diffParticlePosition(0), diffParticlePosition(1), diffParticlePosition(2));
// 			//now, make sure particle position is still on same spherical shell as before
// 			particle_horizontal_movement(*mpas1, particlePosition, diffParticlePosition);

// 			//now can do any vertical movements independent of the horizontal movement
// 			zLevelParticle = zLevelParticle + diffParticlePositionVert;

// 			// dprint(" timestep03 %d, particlePosition (%.7f %.7f %.7f), zLevelParticle %.7f", timeStep, particlePosition[0], particlePosition[1],particlePosition[2], zLevelParticle);



// 			// update segment
// 			Pt p;
// 			p.coords[0] = particlePosition[0]; 
// 			p.coords[1] = particlePosition[1]; 
// 			p.coords[2] = particlePosition[2]; 
// 			seg.pts.push_back(p);

// 		} // timestep loop ends



// 		// push segment to block segment vector
// 		mpas1->segments.push_back(seg);

// 		// updating the mpas1->particles vector for next round
// 		EndPt pt;
// 		pt.pid = mpas1->particles[pid].pid; // Needs modification of diy code to be effective
// 		pt[0] = particlePosition[0];  pt[1] = particlePosition[1];  pt[2] = particlePosition[2];
// 		pt.zLevelParticle = zLevelParticle;
// 		// mpas1->particles[pid].glCellIdx = iCell; 
// 		pt.glCellIdx = iCell; 
// 		particles_finished.push_back(pt);
		

// 		// dprint("particle finished in epoch");
// 	} // particle loop

// 	// dprint("pldone");

// 	mpas1->particles = std::move(particles_finished);

}

void pathline::get_validated_cell_id(const block &mpas1, Array3d &xSubStep, int &iCell, int &nCellVertices)
{

		get_nearby_cell_index_sl(mpas1.nCells, &mpas1.xCell[0], &mpas1.yCell[0], &mpas1.zCell[0], xSubStep(0), xSubStep(1), xSubStep(2), mpas1, iCell, &mpas1.cellsOnCell[0], &mpas1.nEdgesOnCell[0]);

		// // dprint("iCell guess %d", iCell);
		// get_nearby_cell_index(mpas1.nCells, &mpas1.xCell[0], &mpas1.yCell[0], &mpas1.zCell[0], xSubStep(0), xSubStep(1), xSubStep(2), mpas1, iCell, &mpas1.cellsOnCell[0], &mpas1.nEdgesOnCell[0]);

	// if (iCell==-1){
	/* Get cell id */

		// Eigen::VectorXd q(3);
		// q<<xSubStep[0], xSubStep[1], xSubStep[2];//-4677490.455831 4114410.084085 -1336140.586176
		// // q<<-4677490.455831, 4114410.084085, -1336140.586176;
		// // get nearest cell neighbor ids using mpas_c
		// Eigen::VectorXi nearest_cell_idx(1);
		// Eigen::VectorXd dists2_cell(1);
		// mpas1.nns_cells->knn(q, nearest_cell_idx,dists2_cell, 1,0, Nabo::NNSearchF::SORT_RESULTS| Nabo::NNSearchF::ALLOW_SELF_MATCH);
		// iCell = nearest_cell_idx[0] + 1; // cell fortran global id

	// } else{


	// }

	// dprint("iCell found %d", iCell);

	nCellVertices = mpas1.nEdgesOnCell[iCell];
	// dprint("nCellVertices %d", nCellVertices);
}

void pathline::velocity_time_interpolation(const int timeInterpOrder,
										   const double *timeCoeff,
										   const int iCell,
										   const int iLevel,
										   const int nCellVertices,
										   const double zSubStep,
										   block *b,
										   mpas_io &mpas1, // includes verticesOnCell, boundaryVertex, zMid, zTop, xVertex, yVertex, zVertex
										   const Eigen::Array3d &xSubStep,
										   Eigen::Array3d &particleVelocity,
										   double &particleVelocityVert)
{

	double vertCoords[nCellVertices][3];
	double uvCell[nCellVertices][3];
	double verticalVelocityInterp;
	double areaB[nCellVertices];

	

	// get horizontal vertex locations

	for (int aVertex = 0; aVertex < nCellVertices; aVertex++)
	{
		// vertCoords[aVertex][0] = mpas1.xVertex[mpas1.vertexIndex[mpas1.verticesOnCell[(iCell)*mpas1.maxEdges + aVertex]]];
		// vertCoords[aVertex][1] = mpas1.yVertex[mpas1.vertexIndex[mpas1.verticesOnCell[(iCell)*mpas1.maxEdges + aVertex]]];
		// vertCoords[aVertex][2] = mpas1.zVertex[mpas1.vertexIndex[mpas1.verticesOnCell[(iCell)*mpas1.maxEdges + aVertex]]];

		vertCoords[aVertex][0] = b->xVertex[-1+b->verticesOnCell[(iCell)*b->maxEdges + aVertex]];
		vertCoords[aVertex][1] = b->yVertex[-1+b->verticesOnCell[(iCell)*b->maxEdges + aVertex]];
		vertCoords[aVertex][2] = b->zVertex[-1+b->verticesOnCell[(iCell)*b->maxEdges + aVertex]];

		// if (iCell == 5470)
		// 	dprint("vertex on 5470 %d", b->verticesOnCell[(iCell)*b->maxEdges + aVertex]);

	}


	// exit(0);

	// initialize velocities to 0
	particleVelocity << 0,0,0;
	particleVelocityVert = 0;

	// general interpolation for the velocity field
	for (int aTimeLevel = 0; aTimeLevel < timeInterpOrder; aTimeLevel++)
	{

		// get uVertexVelocity, vVertexVelocity, wVertexVelocity for current time: done implicitly by update_velocity_vectors()

		// if (iCell==1000)
		// 	dprint("VONCELL %d %d %d", b->verticesOnCell[iCell * b->maxEdges], b->verticesOnCell[iCell * b->maxEdges +1], b->verticesOnCell[iCell * b->maxEdges +2]);

		// dprint("testval %f", mpas1.zMid[iCell * b->nVertLevels + iLevel]);
		particle_vertical_treatment(nCellVertices, &b->verticesOnCell[iCell * b->maxEdges], &b->uVertexVelocity[0], &b->vVertexVelocity[0], &b->wVertexVelocity[0], b, mpas1, uvCell, iLevel, zSubStep, &b->zMid[iCell * b->nVertLevels], &b->zTop[iCell * b->nVertLevels], &b->vertVelocityTop[iCell * b->nVertLevelsP1], verticalVelocityInterp);

		// dprint("b->uVertexVelocity size %ld, uVertVels %f %f %f %f %f %f", b->uVertexVelocity.size(), b->uVertexVelocity[b->nVertLevels*4894], 
		// b->uVertexVelocity[b->nVertLevels*14825], 
		// b->uVertexVelocity[b->nVertLevels*7728], 
		// b->uVertexVelocity[b->nVertLevels*7727], 
		// b->uVertexVelocity[b->nVertLevels*13860], 
		// b->uVertexVelocity[b->nVertLevels*4895]);

		particleVelocityVert = particleVelocityVert + timeCoeff[aTimeLevel] * verticalVelocityInterp;

		particleVelocity = particleVelocity + timeCoeff[aTimeLevel] * particle_horizontal_interpolation(mpas1, nCellVertices, vertCoords, xSubStep, uvCell, areaB);
	}

	// if (b->verticesOnCell[(iCell)*b->maxEdges + 0] - 1==6464){
	// 	dprint("iCell %d | vel %f %f %f | %f", iCell, particleVelocity(0), particleVelocity(1), particleVelocity(2), particleVelocityVert);
	// }


	// dprint("Exiting"); exit(0);
}

Array3d pathline::particle_horizontal_interpolation(mpas_io &mpas1, const int nCellVertices, const double vertCoords[][3], const Array3d &pointInterp, const double uvCell[][3], double *areaB)
{

	double lambda[nCellVertices];
	wachspress_coordinates(mpas1, nCellVertices, vertCoords, pointInterp, areaB, &lambda[0]);

	// dprint("lambda (%f %f %f %f %f %f) vertCoords %f %f %f", lambda[0], lambda[1], lambda[2],lambda[3], lambda[4], lambda[5], vertCoords[0][0], vertCoords[0][1], vertCoords[0][2]);

	// update particle velocities via horizontal interpolation
	Array3d particle_horizontal_interpolation;
	particle_horizontal_interpolation(0) = wachspress_interpolate(lambda, uvCell, 0, nCellVertices);
	particle_horizontal_interpolation(1) = wachspress_interpolate(lambda, uvCell, 1, nCellVertices);
	particle_horizontal_interpolation(2) = wachspress_interpolate(lambda, uvCell, 2, nCellVertices);

	// dprint("particle_horizontal_interpolation (%f %f %f)", particle_horizontal_interpolation(0), particle_horizontal_interpolation(1), particle_horizontal_interpolation(2));

	return particle_horizontal_interpolation;
}

double pathline::wachspress_interpolate(const double *lambda, const double phi[][3], const int idx, const int nCellVertices)
{

	double sum = 0;

	for (int i = 0; i < nCellVertices; i++)
		sum += lambda[i] * phi[i][idx];

	return sum;
}

void pathline::wachspress_coordinates(const mpas_io &mpas1, const int nVertices, const double vertCoords[][3], const Array3d &pointInterp, double *areaB, double *wachspress, int flag)
{	
	if (flag==1)
	dprint("in wc");
	double areaA[nVertices], wach[nVertices];

	double p[3] = {pointInterp(0), pointInterp(1), pointInterp(2)};
	int im1, i0, ip1;
	for (int i = 0; i < nVertices; i++)
	{
		im1 = (nVertices + i - 2) % nVertices;
		i0 = (nVertices + i - 1) % nVertices;
		ip1 = (nVertices + i) % nVertices;

		double a[3] = {vertCoords[im1][0], vertCoords[im1][1], vertCoords[im1][2]};
		double b[3] = {vertCoords[i0][0], vertCoords[i0][1], vertCoords[i0][2]};
		double c[3] = {vertCoords[ip1][0], vertCoords[ip1][1], vertCoords[ip1][2]};

		areaB[i] = triangle_signed_area(a, b, c, mpas1.radius, flag);
		areaA[i] = triangle_signed_area(p, b, c, mpas1.radius, flag);
		// if (flag==1)
		// 	dprint("areaA[i] %.20g areaB[i] %.20g, a (%f %f %f), b (%f %f %f) c(%f %f %f) p(%f %f %f), im1 %d i0 %d ip1 %d", areaA[i], areaB[i], a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2], p[0], p[1], p[2], im1, i0, ip1);
		// dprint("area %.20g", areaB[i]);
	}

	/* For each vertex compute wachspress coordinate */
	for (int i = 0; i < nVertices; i++)
	{
		wach[i] = areaB[i];
		for (int j = i + 1; j < i + nVertices - 1; j++)
		{

			i0 = (nVertices + j) % nVertices;
			// dprint("%d %d", i, i0);
			//accumulate products for A_ij subareas
			wach[i] = wach[i] * areaA[i0];
		}
		// dprint(" ");
	}

	double wach_total = 0;
	for (int i = 0; i < nVertices; i++)
	{
		wach_total += wach[i];
	}
	// dprint("wach_total %20f", wach_total);

	/* Compute lambda */

	// std::vector<double> wachspress(nVertices);
	wachspress[nVertices - 1] = wach[0] / wach_total;
	for (int i = 1; i < nVertices; i++)
	{
		wachspress[i - 1] = wach[i] / wach_total;
	}
}

void pathline::get_bounding_indices(int &iHigh, int &iLow, const double phiInterp, const double *phiVals, const int iLevel, const int nVertLevels)
{

	if (phiInterp > phiVals[iLevel])
	{
		iHigh = iLevel + 1;
		iLow = iLevel;

		// todo:001 should not come here!!!

		// dprint("Shouldn't be here");
		// assert(0);
	}
	else
	{

		iHigh = iLevel;
		iLow = iLevel + 1;
	}

	// check to make sure points are in range
	// optimization point: smarter algorithm won't have to call this ever
	if (!((phiInterp <= phiVals[iHigh]) && (phiInterp >= phiVals[iLow])))
	{
		// todo:002
		// assert(0);
		get_bounding_indices_brute_force(nVertLevels, phiInterp, phiVals, iLow, iHigh);
	}
}

void pathline::get_bounding_indices_brute_force(const int nVertLevels, const double phiInterp, const double *phiVals, int &iHigh, int &iLow)
{

	// make no assumptions
	for (int aLevel = 0; aLevel < nVertLevels - 1; ++aLevel)
	{
		if (phiVals[aLevel] <= phiInterp && phiInterp <= phiVals[aLevel + 1])
		{
			iLow = aLevel;
			iHigh = aLevel + 1;
			break;
		}
		else if (phiVals[aLevel + 1] <= phiInterp && phiInterp <= phiVals[aLevel])
		{
			iLow = aLevel + 1;
			iHigh = aLevel;
			break;
		}
	}
	// assert(phiVals[iLow] <= phiInterp);
	// assert(phiInterp <= phiVals[iHigh]);
}

void pathline::particle_vertical_treatment(const int nCellVertices, const int *verticesOnCell, const double *uVertexVelocity, const double *vVertexVelocity, const double *wVertexVelocity, block *b, mpas_io &mpas1, double uvCell[][3], const int iLevel, const double zLoc, const double *zMid, const double *zTop, const double *zVertVelocityTop, double &verticalVelocityInterp)
{

	// vertical treatment = passiveFloat

	// dprint("here1 zMid[iLevel] %f", zMid[iLevel]);
	// interpolate vertical velocity
	verticalVelocityInterp = interp_vert_velocity_to_zlevel(b, mpas1, iLevel, zLoc, zTop, zVertVelocityTop);

	// dprint("here2");

	// interpolate the horizontal velocity based on z-levels
	interp_nodal_vectors(b, mpas1, nCellVertices, verticesOnCell, iLevel, b->nVertLevels, zLoc, zMid, uVertexVelocity, vVertexVelocity, wVertexVelocity, uvCell);

	// dprint("here3");

	// dprint("uvCell %f %f %f, %f %f %f", uvCell[0][0], uvCell[0][1], uvCell[0][2], uvCell[1][0], uvCell[1][1], uvCell[1][2]);

	if (verticesOnCell[0] - 1==6464){
		// dprint("vertVelInterp %f iLevel %d vertVelocityTop[iLevel] %f", verticalVelocityInterp, iLevel, zVertVelocityTop[iLevel]);
	}
}

void pathline::interp_nodal_vectors(block *b, mpas_io &mpas1, const int nCellVertices, const int *verticesOnCell, const int iLevel, int nVertLevels, const double phiInterp, const double *phiVals, const double *uVertexVelocity, const double *vVertexVelocity, const double *wVertexVelocity, double uvCell[][3])
{	
	// dprint("uVertexVelocity[0] %f", uVertexVelocity[0]);
	// dprint("iLevell %d, nVertLevels %d", iLevel, nVertLevels);
	int theVertex;

	if (iLevel < 1)
	{
		if (iLevel == 0)
		{
			for (int aVertex = 0; aVertex < nCellVertices; aVertex++)
			{	
				theVertex = verticesOnCell[aVertex] - 1;
				// theVertex = mpas1.vertexIndex[verticesOnCell[aVertex]];
				uvCell[aVertex][0] = uVertexVelocity[theVertex * nVertLevels + 0]; // 0 = maxloc(phiVals(1:nVertLevels),1)
				uvCell[aVertex][1] = vVertexVelocity[theVertex * nVertLevels + 0]; // 0 = maxloc(phiVals(1:nVertLevels),1)
				uvCell[aVertex][2] = wVertexVelocity[theVertex * nVertLevels + 0]; // 0 = maxloc(phiVals(1:nVertLevels),1)
			}
		}
		else if (iLevel == -1)
		{
			for (int aVertex = 0; aVertex < nCellVertices; aVertex++)
			{	
				theVertex = verticesOnCell[aVertex] - 1;
				// theVertex = mpas1.vertexIndex[verticesOnCell[aVertex]];
				uvCell[aVertex][0] = uVertexVelocity[theVertex * nVertLevels + nVertLevels - 1]; // nVertLevels-1 = minloc(..)
				uvCell[aVertex][1] = vVertexVelocity[theVertex * nVertLevels + nVertLevels - 1]; // nVertLevels-1 = minloc(..)
				uvCell[aVertex][2] = wVertexVelocity[theVertex * nVertLevels + nVertLevels - 1]; // nVertLevels-1 = minloc(..)
			}
		}
	}

	int iHigh, iLow;
	const double eps = 1e-14;
	double alpha = 0;
	

	// dprint("iLevel %d, iHigh %d, iLow %d, phiInterp %f", iLevel, iHigh, iLow, phiInterp);
	get_bounding_indices(iHigh, iLow, phiInterp, phiVals, iLevel, nVertLevels);

	// dprint("iLevel %d, iHigh %d, iLow %d, phiInterp %f, phiVals[iLo iHi] %f %f", iLevel, iHigh, iLow, phiInterp, phiVals[iLow], phiVals[iHigh]);
	// dprint("uVertexVelocity[0,1,2] %f %f %f", uVertexVelocity[theVertex*nVertLevels], uVertexVelocity[theVertex*nVertLevels+1], uVertexVelocity[theVertex*nVertLevels+2]);

	// dprint("iLevel %d, iHigh %d, iLow %d, phiInterp %f, phiVals[iLevel] %f", iLevel, iHigh, iLow, phiInterp, phiVals[iLevel]);

	// interpolate to vertical level now
	if (abs(phiVals[iHigh] - phiVals[iLow]) < eps)
	{
		alpha = -0.5;
	}
	else
	{
		alpha = (phiInterp - phiVals[iLow]) / (phiVals[iHigh] - phiVals[iLow]);
	}

	// theVertex = verticesOnCell[0] - 1;
	// int theVertex2 = verticesOnCell[1] - 1;
	// int theVertex3 = verticesOnCell[2] - 1;
	// if (theVertex==6464)
	// 	dprint("theVertex %d, iLow %d, uVertexVelocity[theVertex*nVertLevels + iHigh] %f %f %f, theVertex*nVertLevels + iLow %d", theVertex, iLow, uVertexVelocity[theVertex*nVertLevels + iLow], uVertexVelocity[theVertex2*nVertLevels + iLow], uVertexVelocity[theVertex3*nVertLevels + iLow], (theVertex*nVertLevels + iLow));
	
	// interpolate to the vertical level
	int flag = 0;
	for (int aVertex = 0; aVertex < nCellVertices; ++aVertex)
	{
		theVertex = verticesOnCell[aVertex] - 1;
		// theVertex = mpas1.vertexIndex[verticesOnCell[aVertex]];
		uvCell[aVertex][0] = alpha * uVertexVelocity[theVertex * nVertLevels + iHigh] + (1 - alpha) * uVertexVelocity[theVertex * nVertLevels + iLow];
		uvCell[aVertex][1] = alpha * vVertexVelocity[theVertex * nVertLevels + iHigh] + (1 - alpha) * vVertexVelocity[theVertex * nVertLevels + iLow];
		uvCell[aVertex][2] = alpha * wVertexVelocity[theVertex * nVertLevels + iHigh] + (1 - alpha) * wVertexVelocity[theVertex * nVertLevels + iLow];

		int theVertex2 = verticesOnCell[1] - 1;
		int theVertex3 = verticesOnCell[2] - 1;	
		int theVertex4 = verticesOnCell[3] - 1;	
		int theVertex5 = verticesOnCell[4] - 1;	
		int theVertex6 = verticesOnCell[5] - 1;	
		if (theVertex==6464){
			// dprint(" uVertexVel (%f %f), (%f %f), (%f %f), (%f %f), (%f %f), (%f %f)",  
			// uVertexVelocity[theVertex * nVertLevels + iHigh], uVertexVelocity[theVertex * nVertLevels + iLow],
			// uVertexVelocity[theVertex2 * nVertLevels + iHigh], uVertexVelocity[theVertex2 * nVertLevels + iLow],
			// uVertexVelocity[theVertex3 * nVertLevels + iHigh], uVertexVelocity[theVertex3 * nVertLevels + iLow],
			// uVertexVelocity[theVertex4 * nVertLevels + iHigh], uVertexVelocity[theVertex4 * nVertLevels + iLow],
			// uVertexVelocity[theVertex5 * nVertLevels + iHigh], uVertexVelocity[theVertex5 * nVertLevels + iLow],
			// uVertexVelocity[theVertex6 * nVertLevels + iHigh], uVertexVelocity[theVertex6 * nVertLevels + iLow]
			// );
			flag = 1;
		}

		
	}
	// if (verticesOnCell[0] - 1==6464)
	// 	dprint("uvCell %f %f %f %f %f %f, IH %d IL %d", uvCell[0][0], uvCell[1][0], uvCell[2][0], uvCell[3][0], uvCell[4][0], uvCell[5][0], iHigh, iLow);

	// theVertex = verticesOnCell[0] - 1;
	// int theVertex2 = verticesOnCell[1] - 1;
	// int theVertex3 = verticesOnCell[2] - 1;
	// if (theVertex==6464)
	// 	dprint("uvCell %f %f %f", uvCell[0][0], uvCell[1][0], uvCell[2][0]);
}

double pathline::interp_vert_velocity_to_zlevel(block *b, mpas_io &mpas1, const int iLevel, const double zSubStep, const double *zTop, const double *vertVelocityTop)
{

	double alpha = 0, interp_vert_velocity_to_zlevel;
	// dprint("iLevel %d, alpha %f, vertVelocityTop[iLevel] %f", iLevel, alpha, vertVelocityTop[0]);

	if (iLevel < 1)
	{
		if (iLevel == 0)
		{
			interp_vert_velocity_to_zlevel = vertVelocityTop[0];
		}
		else if (iLevel == -1)
		{
			interp_vert_velocity_to_zlevel = vertVelocityTop[b->nVertLevelsP1 - 1];
		}
	}
	else
	{

		// interpolate the vertical velocity [assumes zTop(iLevel+1) <= zSubStep <= zTop(iLevel)]
		if (zTop[iLevel + 1] <= zSubStep && zSubStep <= zTop[iLevel])
		{
			alpha = (zSubStep - zTop[iLevel + 1]) / (zTop[iLevel] - zTop[iLevel + 1]);
		}
		else if (zTop[iLevel] <= zSubStep && zSubStep <= zTop[iLevel + 1])
		{
			alpha = (zSubStep - zTop[iLevel]) / (zTop[iLevel + 1] - zTop[iLevel]);
		}

		interp_vert_velocity_to_zlevel = alpha * vertVelocityTop[iLevel] + (1 - alpha) * vertVelocityTop[iLevel + 1];
	}
	// dprint("vertVelocityTop iLevel %d %.10f %.10f, %.13f", iLevel, vertVelocityTop[iLevel], vertVelocityTop[iLevel+1], interp_vert_velocity_to_zlevel);
	// dprint("interp_vert_velocity_to_zlevel %f", interp_vert_velocity_to_zlevel);
	return interp_vert_velocity_to_zlevel;
}

void pathline::update_velocity_vectors(mpas_io &mpas1, int frame_no)
{

	// uVertexVelocities[0] = &mpas1.velocityXv[(frame_no + 1) % 2][0];
	// vVertexVelocities[0] = &mpas1.velocityYv[(frame_no + 1) % 2][0];
	// wVertexVelocities[0] = &mpas1.velocityZv[(frame_no + 1) % 2][0];

	// uVertexVelocities[1] = &mpas1.velocityXv[frame_no % 2][0];
	// vVertexVelocities[1] = &mpas1.velocityYv[frame_no % 2][0];
	// wVertexVelocities[1] = &mpas1.velocityZv[frame_no % 2][0];

	// zMid_cur = &mpas1.zMid[(frame_no + 1) % 2][0];
	// zTop_cur = &mpas1.zTop[(frame_no + 1) % 2][0];

	// zMid_cur_size = mpas1.zMid[(frame_no + 1) % 2].size();
	// zTop_cur_size = mpas1.zTop[(frame_no + 1) % 2].size();

	// vertVelocityTop_cur = &mpas1.vertVelocityTop[(frame_no + 1) % 2][0];
	// // dprint("vertVelocityTop_cur size %ld", mpas1.vertVelocityTop[(frame_no + 1) % 2].size());

	// // initialize first frame to zero if first round of advection
	// temp.resize(mpas1.velocityXv[frame_no % 2].size());
	// if (frame_no == 2){
	// 	uVertexVelocities[0] = &temp[0];
	// 	vVertexVelocities[0] = &temp[0];
	// 	wVertexVelocities[0] = &temp[0];
	// }

	// // if (frame_no==2){
	// // 	// uVertexVelocities[0] = &mpas1.velocityXv[(frame_no) % 2][0];
	// // 	// vVertexVelocities[0] = &mpas1.velocityYv[(frame_no) % 2][0];
	// // 	// wVertexVelocities[0] = &mpas1.velocityZv[(frame_no) % 2][0];
	// // 	// zMid_cur = &mpas1.zMid[(frame_no) % 2][0];
	// // 	// zTop_cur = &mpas1.zTop[(frame_no) % 2][0];
	// // 	// zMid_cur_size = mpas1.zMid[(frame_no) % 2].size();
	// // 	// zTop_cur_size = mpas1.zTop[(frame_no) % 2].size();

	// // 	// std::vector<double> temp(mpas1.vertVelocityTop.size());

	// // 	// vertVelocityTop_cur = &mpas1.vertVelocityTop[(frame_no) % 2][0];
	// // 	uVertexVelocities[1] = uVertexVelocities[0]; // &mpas1.velocityXv[(frame_no) % 2][0];
	// // 	vVertexVelocities[1] = vVertexVelocities[0]; // &mpas1.velocityYv[(frame_no) % 2][0];
	// // 	wVertexVelocities[1] = wVertexVelocities[0]; // &mpas1.velocityZv[(frame_no) % 2][0];



	// // }
}
