#include "pathline.h"
#include "mpas_io.h"
#include "misc.h"
#include <Eigen/Dense>
#include "geometry_utils.h"
#include "block.h"
#include "advect.h"
#include <cmath>

using namespace Eigen;

bool pathline::in_global_domain(const Pt& p){

	if (pow(-6352390-p.coords[0],2)+pow(-447879-p.coords[1],2)+(pow(197714-p.coords[2],2))>pow(5524190.98744,2))
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
	dprint("bgid %d iCell %d b->in_partition[iCell] %d round %d  b->gcIdxToGid[iCell] %d", b->gid, iCell, b->in_partition[iCell], round,  b->gcIdxToGid[iCell]);
	if (b->in_partition[iCell] == round)
		return true;
	else
		return false;
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

bool pathline::compute_streamlines(block *b, const diy::Master::ProxyWithLink &cp, const diy::Assigner &assigner){

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

	// dprint("mpas1->particles.size() %ld, nCells %ld", mpas1->particles.size(), mpas1->nCells);

	
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

						// get new time step (tm = (timestep-1)*dt)
						tSubStep = (timeStep + kWeightT(subStep)) * dt;

						get_validated_cell_id(*b, xSubStep, iCell, nCellVertices);
						iLevel = get_vertical_id(b->maxLevelCell[iCell], zSubStep, &zMid_cur[iCell * b->nVertLevels]);
						dprint("iCell %d, iLevel %d, maxLevelCell %d", iCell, iLevel, b->maxLevelCell[iCell]);

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
											*b,
											xSubStep,
											particleVelocity,
											particleVelocityVert);
						dprint("vel %d %d| %f %f %f", iCell, iLevel, particleVelocity(0), particleVelocity(1), particleVelocity(2));
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
						seg.pts.push_back(p);

						cur_nsteps ++;

						// check if particle inside global domain and handle
						if (!in_global_domain(p) || cur_nsteps >= nSteps){
							finished = true;
							dprint("particle exited gid %d pid %d, cur_nsteps %d, in_global_domain(p) %d, nSteps %f", b->gid, b->particles[pid].pid, cur_nsteps, in_global_domain(p), nSteps);
							break;
						}

						// check if inside local domain and handle
						int round = 1;
						if (!in_local_domain(b , p, iCell, round)){

							dprint("jumped!! in %d, cur_nsteps %d, pid %d", b->gid, cur_nsteps, b->particles[pid].pid);
							break;
						}



			} // timeStep loop ends


			// if unfinished, enqueue 
			if (finished == false){
				int dest_gid = b->gcIdxToGid[iCell];
				int dest_proc = assigner.rank(dest_gid);
				diy::BlockID dest_block = {dest_gid, dest_proc};
				dprint("iCell %d, sou %d des %d %d, cur_nsteps %d", iCell, b->gid, dest_gid, dest_proc, cur_nsteps);
				EndPt pt;
				pt.pid = b->particles[pid].pid; // Needs modification of diy code to be effective
				pt[0] = particlePosition[0];  pt[1] = particlePosition[1];  pt[2] = particlePosition[2];
				pt.nsteps = b->particles[pid].nsteps;
				pt.zLevelParticle = zLevelParticle;
				// mpas1->particles[pid].glCellIdx = iCell; 
				pt.glCellIdx = iCell;
				cp.enqueue(dest_block, pt);
			}

			// push segment to block segment vector
			// dprint("segsize %ld", seg.pts.size());
			b->segments.push_back(seg);
			dprint("numsegs %ld", b->segments.size());

		} // particle loop

	
	
	
	if (b->particles_store.size() > 0 ){
		dprint("moving from store in %d", b->gid);
		b->particles = std::move(b->particles_store);
		// dprint("finishedt a callback in gid %d, %ld %ld", b->gid,  b->particles.size(), b->particles_store.size());
		return false;
	}
	else {
		// b->particle_store.clear();
		dprint("finishedf a callback in gid %d, %ld %ld", b->gid,  b->particles.size(), b->particles_store.size());
		b->particles.clear(); // clearing particles if no particles came in
		return true;
	}
		
	

}

void pathline::compute_epoch(block *mpas1, int framenum)
{

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

	size_t nCellsnVertLeves = mpas1->nCells * mpas1->nVertLevels;

	// dprint("mpas1->particles.size() %ld, nCells %ld", mpas1->particles.size(), mpas1->nCells);

	std::vector<EndPt> particles_finished;

	for (int pid = 0; pid < mpas1->particles.size(); pid++)
	{

		// dprint("pid %d", pid);
		// create segment
		Segment seg;
		seg.pid = mpas1->particles[pid].pid;

		// if (pid != 8 ) continue;
		// if (mpas1->gid != 1 ) continue;

		// Array3d particlePosition = {
		// 	particles_current[3*pid+0],
		// 	particles_current[3*pid+1],
		// 	particles_current[3*pid+2]};

		Array3d particlePosition = {
			mpas1->particles[pid][0],
			mpas1->particles[pid][1],
			mpas1->particles[pid][2]};

		// read zLevelParticle
		// zLevelParticle = mpas1->zLevelParticle[simStep*nParticles+pid];
		// zLevelParticle = zparticle_current[pid];
		zLevelParticle = mpas1->particles[pid].zLevelParticle; // particles_zLevel[0];
		// iCell = -1;											   //mpas1->particles[pid].glCellIdx;

		

		iCell = mpas1->particles[pid].glCellIdx;

		// dprint("pid %d zLevelParticle %f, particlePosition (%f %f %f)", pid, zLevelParticle, particlePosition[0], particlePosition[1], particlePosition[2]);

		// dprint("zLevelParticle %f", zLevelParticle);
		bool finished = false;
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
					particle_horizontal_movement(*mpas1, xSubStep, diffSubStep);
				}

				// vertical
				zSubStep = zLevelParticle;
				diffSubStepVert = kWeightKVert(subStep) * kCoeffVert(subStep);

				// dprint("timestep00 %d, subStep %d, pid %d, (%.7f %.7f %.7f), diff % f (%f %f %f) zLevel %.7f", timeStep, subStep, mpas1->particles[pid].pid, xSubStep(0), xSubStep(1), xSubStep(2), diffSubStepVert, diffSubStep(0), diffSubStep(1), diffSubStep(2), zLevelParticle);


				if (kWeightKVert(subStep) != 0)
				{
					zSubStep = zSubStep + diffSubStepVert;
				}

				// get new time step (tm = (timestep-1)*dt)
				tSubStep = (timeStep + kWeightT(subStep)) * dt;

				get_validated_cell_id(*mpas1, xSubStep, iCell, nCellVertices);

				// dprint("pid %d, Cell %d, gid %d", pid, iCell, mpas1->gid);
				

				// get cellIdx

				iCell = mpas1->cellIndex[iCell];
				

				// dprint("iCell %d %f", iCell, mpas1->xCell[iCell]);

				// dprint("zMid_cur %ld, zTop_cur %ld", zMid_cur_size, zTop_cur_size);
				iLevel = get_vertical_id(mpas1->maxLevelCell[iCell], zSubStep, &zMid_cur[iCell * mpas1->nVertLevels]);
				// dprint("iLevel %d", iLevel);

				int timeInterpOrder = 2;
				double timeCoeff[2];
				timeCoeff[0] = tSubStep / dtSim;
				timeCoeff[1] = 1.0 - timeCoeff[0];
				// dprint("mid");
				// dprint("tSubStep %f, timeCoeff (%f %f)", tSubStep, timeCoeff[0], timeCoeff[1]);

				// velocity_time_interpolation
				// double v[3];
				// double p[3] = {xSubStep(0), xSubStep(1), xSubStep(2)};

				velocity_time_interpolation(timeInterpOrder,
											timeCoeff,
											iCell,
											iLevel,
											nCellVertices,
											zSubStep,
											*mpas1,
											xSubStep,
											particleVelocity,
											particleVelocityVert);

				// dprint(" %f (%f %f %f)", zLevelParticle, particleVelocity[0], particleVelocity[1],particleVelocity[2]); 
				// dprint("timestep01 %d,subStep %d, pid %d,(%f %f %f), velVert %.07e, particleVelocity (%.07e %.7ef %.7e), zLevel %f", timeStep, subStep, mpas1->particles[pid].pid, xSubStep(0), xSubStep(1), xSubStep(2), particleVelocityVert, particleVelocity(0), particleVelocity(1), particleVelocity(2), zLevelParticle);


				//!!!!!!!!!! FORM INTEGRATION WEIGHTS kj !!!!!!!!!!
				kCoeff.col(subStep + 1) = dt * particleVelocity;
				kCoeffVert(subStep + 1) = dt * particleVelocityVert;

				// dprint(" %f (%f )", zLevelParticle, kCoeff.col(subStep + 1)(0)); 
				iCell = mpas1->indexToCellID[iCell];
			} // substep loop ends

			

			// update particle positions
			diffParticlePosition << 0,0,0;
			diffParticlePositionVert = 0;

			for (int subStep = 0; subStep < subStepOrder; subStep++){

              		//first complete particle integration
					diffParticlePosition = diffParticlePosition + kWeightX.col(subStep) * kCoeff.col(subStep+1);
					diffParticlePositionVert = diffParticlePositionVert + kWeightXVert(subStep) * kCoeffVert(subStep+1);
			}

			// dprint("timestep02 diffParticlePosition %.7e %.7e %.7e", diffParticlePosition(0), diffParticlePosition(1), diffParticlePosition(2));
			//now, make sure particle position is still on same spherical shell as before
			particle_horizontal_movement(*mpas1, particlePosition, diffParticlePosition);

			//now can do any vertical movements independent of the horizontal movement
			zLevelParticle = zLevelParticle + diffParticlePositionVert;

			// dprint(" timestep03 %d, particlePosition (%.7f %.7f %.7f), zLevelParticle %.7f", timeStep, particlePosition[0], particlePosition[1],particlePosition[2], zLevelParticle);



			// update segment
			Pt p;
			p.coords[0] = particlePosition[0]; 
			p.coords[1] = particlePosition[1]; 
			p.coords[2] = particlePosition[2]; 
			seg.pts.push_back(p);

		} // timestep loop ends



		// push segment to block segment vector
		mpas1->segments.push_back(seg);

		// updating the mpas1->particles vector for next round
		EndPt pt;
		pt.pid = mpas1->particles[pid].pid; // Needs modification of diy code to be effective
		pt[0] = particlePosition[0];  pt[1] = particlePosition[1];  pt[2] = particlePosition[2];
		pt.zLevelParticle = zLevelParticle;
		// mpas1->particles[pid].glCellIdx = iCell; 
		pt.glCellIdx = iCell; 
		particles_finished.push_back(pt);
		

		// dprint("particle finished in epoch");
	} // particle loop

	// dprint("pldone");

	mpas1->particles = std::move(particles_finished);

}

void pathline::get_validated_cell_id(const block &mpas1, Array3d &xSubStep, int &iCell, int &nCellVertices)
{
		get_nearby_cell_index_sl(mpas1.nCells, &mpas1.xCell[0], &mpas1.yCell[0], &mpas1.zCell[0], xSubStep(0), xSubStep(1), xSubStep(2), mpas1, iCell, &mpas1.cellsOnCell[0], &mpas1.nEdgesOnCell[0]);

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
}

void pathline::velocity_time_interpolation(const int timeInterpOrder,
										   const double *timeCoeff,
										   const int iCell,
										   const int iLevel,
										   const int nCellVertices,
										   const double zSubStep,
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

		vertCoords[aVertex][0] = mpas1.xVertex[mpas1.verticesOnCell[(iCell)*mpas1.maxEdges + aVertex]-1];
		vertCoords[aVertex][1] = mpas1.yVertex[mpas1.verticesOnCell[(iCell)*mpas1.maxEdges + aVertex]-1];
		vertCoords[aVertex][2] = mpas1.zVertex[mpas1.verticesOnCell[(iCell)*mpas1.maxEdges + aVertex]-1];

		// set areaB[aVertex] here later
		// dprint("iCell %d, vertexIndex[] %d", iCell, mpas1.vertexIndex[mpas1.verticesOnCell[(iCell)*mpas1.maxEdges+aVertex]]);
		// dprint("%f %f %f ", vertCoords[aVertex][0], vertCoords[aVertex][1], vertCoords[aVertex][2]);
	}

	// exit(0);

	// initialize velocities to 0
	particleVelocity << 0,0,0;
	particleVelocityVert = 0;

	// general interpolation for the velocity field

	// for (int aTimeLevel = 0; aTimeLevel < timeInterpOrder; aTimeLevel++)
	// {

		// get uVertexVelocity, vVertexVelocity, wVertexVelocity for current time: done implicitly by update_velocity_vectors()


		// particle_vertical_treatment(nCellVertices, &mpas1.verticesOnCell[iCell * mpas1.maxEdges], uVertexVelocities[aTimeLevel], vVertexVelocities[aTimeLevel], wVertexVelocities[aTimeLevel], mpas1, uvCell, iLevel, zSubStep, &zMid_cur[iCell * mpas1.nVertLevels], &zTop_cur[iCell * mpas1.nVertLevels], &vertVelocityTop_cur[iCell * mpas1.nVertLevelsP1], verticalVelocityInterp);

		// particleVelocityVert = particleVelocityVert + timeCoeff[aTimeLevel] * verticalVelocityInterp;

		// particleVelocity = particleVelocity + timeCoeff[aTimeLevel] * particle_horizontal_interpolation(mpas1, nCellVertices, vertCoords, xSubStep, uvCell, areaB);
		int aTimeLevel = 0;
		particle_vertical_treatment(nCellVertices, &mpas1.verticesOnCell[iCell * mpas1.maxEdges], uVertexVelocities[aTimeLevel], vVertexVelocities[aTimeLevel], wVertexVelocities[aTimeLevel], mpas1, uvCell, iLevel, zSubStep, &zMid_cur[iCell * mpas1.nVertLevels], &zTop_cur[iCell * mpas1.nVertLevels], &vertVelocityTop_cur[iCell * mpas1.nVertLevelsP1], verticalVelocityInterp);

		particleVelocityVert = particleVelocityVert +  verticalVelocityInterp;

		particleVelocity = particleVelocity + particle_horizontal_interpolation(mpas1, nCellVertices, vertCoords, xSubStep, uvCell, areaB);

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

void pathline::wachspress_coordinates(const mpas_io &mpas1, const int nVertices, const double vertCoords[][3], const Array3d &pointInterp, double *areaB, double *wachspress)
{

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

		areaB[i] = triangle_signed_area(a, b, c, mpas1.radius);
		areaA[i] = triangle_signed_area(p, b, c, mpas1.radius);

		// dprint("areaA[i] %.20g areaB[i] %.20g", areaA[i], areaB[i]);
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

void pathline::particle_vertical_treatment(const int nCellVertices, const int *verticesOnCell, const double *uVertexVelocity, const double *vVertexVelocity, const double *wVertexVelocity, mpas_io &mpas1, double uvCell[][3], const int iLevel, const double zLoc, const double *zMid, const double *zTop, const double *zVertVelocityTop, double &verticalVelocityInterp)
{

	// vertical treatment = passiveFloat

	
	// interpolate vertical velocity
	verticalVelocityInterp = interp_vert_velocity_to_zlevel(mpas1, iLevel, zLoc, zTop, zVertVelocityTop);

	// interpolate the horizontal velocity based on z-levels
	interp_nodal_vectors(mpas1, nCellVertices, verticesOnCell, iLevel, mpas1.nVertLevels, zLoc, zMid, uVertexVelocity, vVertexVelocity, wVertexVelocity, uvCell);

	// dprint("uvCell %f %f %f, %f %f %f", uvCell[0][0], uvCell[0][1], uvCell[0][2], uvCell[1][0], uvCell[1][1], uvCell[1][2]);
	
	dprint("verticesOnCell %d %d %d %d %d %d| maxEdges %d| zLoc %f", verticesOnCell[0], verticesOnCell[1], verticesOnCell[2], verticesOnCell[3], verticesOnCell[4], verticesOnCell[5], mpas1.maxEdges, zLoc);
}

void pathline::interp_nodal_vectors(mpas_io &mpas1, const int nCellVertices, const int *verticesOnCell, const int iLevel, int nVertLevels, const double phiInterp, const double *phiVals, const double *uVertexVelocity, const double *vVertexVelocity, const double *wVertexVelocity, double uvCell[][3])
{	

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

	dprint("iLevel %d, iHigh %d, iLow %d, phiInterp %f, phiVals[iLo iHi] %f %f", iLevel, iHigh, iLow, phiInterp, phiVals[iLow], phiVals[iHigh]);
	// dprint("uVertexVelocity[0,1,2] %f %f %f", uVertexVelocity[theVertex*nVertLevels], uVertexVelocity[theVertex*nVertLevels+1], uVertexVelocity[theVertex*nVertLevels+2]);


	// dprint("iLevel %d, iHigh %d, iLow %d, phiInterp %f, phiVals[iLevel] %d", iLevel, iHigh, iLow, phiInterp, 0);

	// interpolate to vertical level now
	if (abs(phiVals[iHigh] - phiVals[iLow]) < eps)
	{
		alpha = -0.5;
	}
	else
	{
		alpha = (phiInterp - phiVals[iLow]) / (phiVals[iHigh] - phiVals[iLow]);
	}

	// dprint("val %f", *uVertexVelocity);
	// interpolate to the vertical level
	for (int aVertex = 0; aVertex < nCellVertices; ++aVertex)
	{
		theVertex = verticesOnCell[aVertex] - 1;
		// theVertex = mpas1.vertexIndex[verticesOnCell[aVertex]];
		uvCell[aVertex][0] = alpha * uVertexVelocity[theVertex * nVertLevels + iHigh] + (1 - alpha) * uVertexVelocity[theVertex * nVertLevels + iLow];
		uvCell[aVertex][1] = alpha * vVertexVelocity[theVertex * nVertLevels + iHigh] + (1 - alpha) * vVertexVelocity[theVertex * nVertLevels + iLow];
		uvCell[aVertex][2] = alpha * wVertexVelocity[theVertex * nVertLevels + iHigh] + (1 - alpha) * wVertexVelocity[theVertex * nVertLevels + iLow];

		dprint("alpha %f iLow %d, iHigh %d, uVertexVelocity[%d]= %f", alpha, iLow, iHigh, theVertex, uVertexVelocity[theVertex * nVertLevels + iLow]);

		// dprint("theVertex %d, iHigh %d, uVertexVelocity[theVertex*nVertLevels + iHigh] %f", theVertex, iHigh, uVertexVelocity[theVertex*nVertLevels + iHigh]);
	}
}

double pathline::interp_vert_velocity_to_zlevel(mpas_io &mpas1, const int iLevel, const double zSubStep, const double *zTop, const double *vertVelocityTop)
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
			interp_vert_velocity_to_zlevel = vertVelocityTop[mpas1.nVertLevelsP1 - 1];
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

	uVertexVelocities[0] = &mpas1.velocityXv[(frame_no + 1) % 2][0];
	vVertexVelocities[0] = &mpas1.velocityYv[(frame_no + 1) % 2][0];
	wVertexVelocities[0] = &mpas1.velocityZv[(frame_no + 1) % 2][0];

	uVertexVelocities[1] = &mpas1.velocityXv[frame_no % 2][0];
	vVertexVelocities[1] = &mpas1.velocityYv[frame_no % 2][0];
	wVertexVelocities[1] = &mpas1.velocityZv[frame_no % 2][0];

	zMid_cur = &mpas1.zMid[(frame_no + 1) % 2][0];
	zTop_cur = &mpas1.zTop[(frame_no + 1) % 2][0];

	zMid_cur_size = mpas1.zMid[(frame_no + 1) % 2].size();
	zTop_cur_size = mpas1.zTop[(frame_no + 1) % 2].size();

	vertVelocityTop_cur = &mpas1.vertVelocityTop[(frame_no + 1) % 2][0];
	// dprint("vertVelocityTop_cur size %ld", mpas1.vertVelocityTop[(frame_no + 1) % 2].size());

	// initialize first frame to zero if first round of advection
	temp.resize(mpas1.velocityXv[frame_no % 2].size());
	if (frame_no == 2){
		uVertexVelocities[0] = &temp[0];
		vVertexVelocities[0] = &temp[0];
		wVertexVelocities[0] = &temp[0];
	}

	// if (frame_no==2){
	// 	// uVertexVelocities[0] = &mpas1.velocityXv[(frame_no) % 2][0];
	// 	// vVertexVelocities[0] = &mpas1.velocityYv[(frame_no) % 2][0];
	// 	// wVertexVelocities[0] = &mpas1.velocityZv[(frame_no) % 2][0];
	// 	// zMid_cur = &mpas1.zMid[(frame_no) % 2][0];
	// 	// zTop_cur = &mpas1.zTop[(frame_no) % 2][0];
	// 	// zMid_cur_size = mpas1.zMid[(frame_no) % 2].size();
	// 	// zTop_cur_size = mpas1.zTop[(frame_no) % 2].size();

	// 	// std::vector<double> temp(mpas1.vertVelocityTop.size());

	// 	// vertVelocityTop_cur = &mpas1.vertVelocityTop[(frame_no) % 2][0];
	// 	uVertexVelocities[1] = uVertexVelocities[0]; // &mpas1.velocityXv[(frame_no) % 2][0];
	// 	vVertexVelocities[1] = vVertexVelocities[0]; // &mpas1.velocityYv[(frame_no) % 2][0];
	// 	wVertexVelocities[1] = wVertexVelocities[0]; // &mpas1.velocityZv[(frame_no) % 2][0];



	// }
}

// setting velocity vectors for streamline computation
void pathline::set_velocity_vectors(mpas_io &mpas1){

	uVertexVelocities[0] = mpas1.uVertexVelocity.data();
	vVertexVelocities[0] = mpas1.vVertexVelocity.data();
	wVertexVelocities[0] = mpas1.wVertexVelocity.data();

	uVertexVelocities[1] = mpas1.uVertexVelocity.data();
	vVertexVelocities[1] = mpas1.vVertexVelocity.data();
	wVertexVelocities[1] = mpas1.wVertexVelocity.data();

	zMid_cur = mpas1.zMid_.data();
	zTop_cur = mpas1.zTop_.data();

	zMid_cur_size = mpas1.zMid_.size();
	zTop_cur_size = mpas1.zTop_.size();

	vertVelocityTop_cur = mpas1.vertVelocityTop_.data();


}