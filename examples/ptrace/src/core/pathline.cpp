#include "pathline.h"
#include "mpas_io.h"
#include "misc.h"
#include <Eigen/Dense>
#include "geometry_utils.h"
#include "block.h"
#include "advect.h"
#include <unordered_set>

using namespace Eigen;


bool pathline::compute_flow(block *b, const diy::Master::ProxyWithLink &cp, const diy::Assigner &assigner, int prediction, size_t &nsteps, std::vector<EndPt> &particles_hold, int skip_rate){


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
	// dprint("b->particles.size() %ld, nCells %ld", b->particles.size(), b->nCells);
	int iCell_prev;

	// std::vector<EndPt> particles_finished;
	// dprint("particles %ld", b->particles.size());

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
			bool finished_final = false;

			
			// if (iCell==1780)
			if (b->particles[pid].pid==1780)
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

						// int timeInterpOrder = 1;
						// double timeCoeff[1];
						// timeCoeff[0] = 1;
						velocity_time_interpolation(timeInterpOrder,
											timeCoeff,
											iCell,
											iLevel,
											nCellVertices,
											zSubStep,
											b,
											xSubStep,
											particleVelocity,
											particleVelocityVert);

						// if (b->particles[pid].pid==3710)
						// 	dprint("vel %d %d| %f %f %f", iCell, iLevel, particleVelocity(0), particleVelocity(1), particleVelocity(2));

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
						// if (prediction == true && cur_nsteps % skip_rate == 0)
						// {   
						//     EndPt way_pt;
						//     way_pt[0] = b->xCell[iCell];
						//     way_pt[1] = b->yCell[iCell];
						//     way_pt[2] = b->zCell[iCell];
						//     way_pt.predonly = 1;
						//     particles_hold.push_back(way_pt);
						// }

						if (iCell_prev == 52571 || iCell_prev == 135779 || iCell_prev == 151446 || iCell_prev == 172732 || iCell_prev == 189476 || iCell_prev == 115789 || iCell_prev == 201643 || iCell_prev == 158103 || iCell_prev == 191753 || iCell_prev == 30044 || iCell_prev == 51769 || iCell_prev == 126417 || iCell_prev == 54758 || iCell_prev == 76847)
						{
							// dprint("particle stopped gid %d, pid %d steps %d", b->gid, b->particles[pid].pid, b->particles[pid].nsteps);
							finished = true;
							break;
						}

						
						// if (iCell==1780)
						if (b->particles[pid].pid==1780)
							dprint("next particle gid %d pid %d, cur_nsteps %d, nSteps %f | iCell %d, (%f %f %f), zLevP %f, iLevel %d", b->gid, b->particles[pid].pid, cur_nsteps, nSteps, iCell, particlePosition(0), particlePosition(1), particlePosition(2), zLevelParticle, iLevel);

						// if (b->particles[pid].pid==3710)
						// 	dprint("next particle gid %d pid %d, cur_nsteps %d, nSteps %f | iCell %d, (%f %f %f), zLevP %f", b->gid, b->particles[pid].pid, cur_nsteps, nSteps, iCell, particlePosition(0), particlePosition(1), particlePosition(2), zLevelParticle);	

						// check if particle inside global domain and handle
						if (!in_global_domain(p, b->xCell[iCell], b->yCell[iCell], b->zCell[iCell])){
						// if (cur_nsteps >= nSteps || iCell == 31823)
						
							finished = true;
							finished_final = true;
							// if (b->particles[pid].pid==5550)
							// 	dprint("exited global %d in gid %d, cur_nsteps %d, nSteps %f", b->particles[pid].pid, b->gid, cur_nsteps, nSteps);
						
							
						}
						
						if (cur_nsteps >= nSteps){
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
						if (!in_local_domain(b , p, iCell, round, b->particles[pid].pid)){
				            if (b->particles[pid].pid==5550)
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
						dprint("alert! jump gid %d dest_gid %d, pid %d, prev iCell %d, iCell %d, cur_nsteps %d, (%f %f %f)", b->gid, dest_gid, b->particles[pid].pid, b->particles[pid].glCellIdx, iCell, cur_nsteps, b->particles[pid][0],  b->particles[pid][1],  b->particles[pid][2]);
					}
			}else if (cur_nsteps == nSteps && finished_final == false){
				// particle continuing to next epoch in the same cell
				EndPt pt;
				pt.pid = b->particles[pid].pid; // Needs modification of diy code to be effective
				pt[0] = particlePosition[0];  pt[1] = particlePosition[1];  pt[2] = particlePosition[2];
				pt.nsteps = 0; // b->particles[pid].nsteps;
				pt.zLevelParticle = zLevelParticle;
				pt.glCellIdx = iCell;

				b->particles_continuing.push_back(pt);
			}

			// push segment to block segment vector
			if (b->particles[pid].pid==1780)
				dprint("segsize %ld, gid %d, segments %ld", seg.pts.size(), b->gid, b->segments.size());
			// dprint ("seg.nsteps %d", seg.nsteps);
			b->segments.push_back(seg);
			// dprint("numsegs %ld", b->segments.size());

		} // particle loop

		// dprint("b->segments size %ld , gid %d", b->segments.size(), cp.gid());
	

	if (b->particles_store.size() > 0 ){
		// dprint("moving from store in %d", b->gid);
		b->particles = std::move(b->particles_store); // incoming particles from different cell in same epoch
		// dprint("finishedt a callback in gid %d, %ld %ld", b->gid,  b->particles.size(), b->particles_store.size());
		return false;
	}
	else {
		b->particles_store.clear();
		// dprint("finishedf a callback in gid %d, %ld %ld", b->gid,  b->particles.size(), b->particles_store.size());
		b->particles.clear(); // clearing particles if no particles came in
		// dprint("b->particles_continuing %ld, gid %d", b->particles_continuing.size(), cp.gid());
		// b->particles = std::move(b->particles_continuing); // to be continued advecting in the next epoch
		
		return true;
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

bool pathline::in_local_domain (const block *b, const Pt& p, int &iCell, int round, int pid){

		int nCellVertices;
		Array3d xSubStep;
		xSubStep(0) = p.coords[0];
		xSubStep(1) = p.coords[1];
		xSubStep(2) = p.coords[2];
		get_validated_cell_id(*b, xSubStep, iCell, nCellVertices);
		
		// if (pid == 5550)
		// 	dprint("bgid %d iCell %d round %d  b->gcIdToGid[iCell] %d", b->gid, iCell, round,  b->gcIdToGid[iCell+1]);
		
		// if (b->gcIdxToGid[iCell] == b->gid)
		// if (b->in_partition[iCell] == round)
		if (b->gcIdToGid[iCell+1] == b->gid)
			return true;
		else
			return false;
}