#include "streamline.h"
#include "block.h"
#include "misc.h"
#include "advect.h"
#include "geometry_utils.h"

using namespace Eigen;

streamline::streamline(mpas_io &mpas1, double dtSim_in, double dtParticle_in)
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


bool streamline::compute_flow(block *b, mpas_io &mpas1, const diy::Master::ProxyWithLink &cp, const diy::Assigner &assigner, int prediction, size_t &nsteps, std::vector<EndPt> &particles_hold, int skip_rate){

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

						// int timeInterpOrder = 2;
						// double timeCoeff[2];
						// timeCoeff[0] = tSubStep / dtSim;
						// timeCoeff[1] = 1.0 - timeCoeff[0];

						int timeInterpOrder = 1;
						double timeCoeff[1];
						timeCoeff[0] = 1;
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