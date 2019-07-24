#include "pathline.h"
#include "mpas_io.h"
#include "misc.h"
#include <Eigen/Dense>
#include "geometry_utils.h"
#include "block.h"
#include "advect.h"

using namespace Eigen;

pathline::pathline(mpas_io &mpas1, double dtSim_in, double dtParticle_in){

	uVertexVelocities.resize(timeInterpOrder);
	vVertexVelocities.resize(timeInterpOrder);
	wVertexVelocities.resize(timeInterpOrder);


	dtSim = dtSim_in;
	dtParticle = dtParticle_in;

	// adjust time step for consistency with integer number of steps
	nSteps = ceil(dtSim/dtParticle);
	dt = dtSim/nSteps;




}





void pathline::compute_epoch(block *mpas1){



	Array2d kWeightK = {0,0.5}, kWeightT = {0, 0.5};
	Array2d kWeightKVert = kWeightK, kWeightTVert = kWeightT;

	ArrayXXd kWeightX(3,2);
	kWeightX << 0, 1,
	0, 1, 
	0, 1;
	int subStepOrder=2;
	Array2d kWeightXVert = kWeightX.block<1,2>(0,0); // see https://eigen.tuxfamily.org/dox/group__TutorialBlockOperations.html
	// Array2d kWeightXVert{0,1};

	// local variables
	Array3d particleVelocity, xSubStep, diffSubStep, diffParticlePosition;
	double tSubStep, diffParticlePositionVert;
	double zSubStep, zLevelParticle, diffSubStepVert;
	int iCell=-1, nCellVertices, iLevel; 
	double particleVelocityVert;

	size_t nCellsnVertLeves = mpas1->nCells*mpas1->nVertLevels;

	dprint("nCellsnVertLeves %ld", nCellsnVertLeves);
	std::vector<EndPt> particles_finished;

	for (int pid=0; pid<mpas1->particles.size(); pid++){
		
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
			zLevelParticle = mpas1->particles[pid].zLevelParticle;// particles_zLevel[0];
			iCell = mpas1->particles[pid].glCellIdx;

			// dprint("simStep %d pid %d zLevelParticle %f, particlePosition (%f %f %f)", simStep, pid, zLevelParticle, particlePosition[0], particlePosition[1], particlePosition[2]);

			// dprint("zLevelParticle %f", zLevelParticle);
			bool finished = false;
			for (int timeStep=0; timeStep<nSteps; timeStep++){

				ArrayXXd kCoeff = ArrayXXd::Zero(3, subStepOrder+1);
				ArrayXXd kCoeffVert = ArrayXXd::Zero(1, subStepOrder+1);

				for (int subStep=0; subStep<subStepOrder; subStep++){

					xSubStep  = particlePosition;
					diffSubStep = kWeightK(subStep)*kCoeff.col(subStep);


					if (kWeightK(subStep) != 0){
						particle_horizontal_movement(*mpas1, xSubStep, diffSubStep);
					}

					

					
					// vertical
					zSubStep = zLevelParticle;
					diffSubStepVert = kWeightKVert(subStep) * kCoeffVert(subStep);

					if (kWeightKVert(subStep) != 0){
						zSubStep = zSubStep + diffSubStepVert;
					}

					// get new time step (tm = (timestep-1)*dt)
					tSubStep = (timeStep + kWeightT(subStep)) * dt;

					get_validated_cell_id(*mpas1, xSubStep, iCell, nCellVertices);

					// iLevel = get_vertical_id(mpas1->maxLevelCell[iCell], zSubStep, &zMid_cur[nCellsnVertLeves+iCell*mpas1->nVertLevels]);


					
					int timeInterpOrder = 2;
					double timeCoeff[2];
					timeCoeff[0] = tSubStep / dtSim;
					timeCoeff[1] = 1.0 - timeCoeff[0];

					// dprint("tSubStep %f, timeCoeff (%f %f)", tSubStep, timeCoeff[0], timeCoeff[1]);

            		// velocity_time_interpolation
					// double v[3]; 
					double p[3] = {xSubStep(0), xSubStep(1), xSubStep(2)};



					//!!!!!!!!!! FORM INTEGRATION WEIGHTS kj !!!!!!!!!!
					kCoeff.col(subStep+1) = dt * particleVelocity;
					kCoeffVert(subStep+1) = dt * particleVelocityVert;
					
				}	


				/* if not in current domain enqueue for exchange break, set flag */
			}

			/* if flag unset, pushback point to particles_finished */

		}

		/* move particles_finished to mpas1->particles */

		/* exchange */

		/* pushback incoming particles to mpas1->particles */

		/* do a reduce of remaining particles */

}


void pathline::get_validated_cell_id(const block &mpas1, Array3d &xSubStep, int &iCell, int &nCellVertices){

				/* Get cell id */
				// Eigen::VectorXd q(3);
				// q<<xSubStep[0], xSubStep[1], xSubStep[2];
				// // get nearest cell neighbor ids using mpas_c
				// Eigen::VectorXi nearest_cell_idx(1);
				// Eigen::VectorXd dists2_cell(1);
				// mpas1.nns_cells->knn(q, nearest_cell_idx,dists2_cell, 1,0, Nabo::NNSearchF::SORT_RESULTS| Nabo::NNSearchF::ALLOW_SELF_MATCH);
				// iCell = nearest_cell_idx[0] ; // cell local id


				get_nearby_cell_index(mpas1.nCells, &mpas1.xCell[0], &mpas1.yCell[0], &mpas1.zCell[0], xSubStep(0), xSubStep(1), xSubStep(2), mpas1, iCell, &mpas1.cellsOnCell[0], &mpas1.nEdgesOnCell[0]);




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
									double &particleVelocityVert){




		double vertCoords[nCellVertices][3];
		double uvCell[nCellVertices][3];
		double verticalVelocityInterp;
		double areaB[nCellVertices];
		

		// get horizontal vertex locations

		for (int aVertex=0; aVertex<nCellVertices; aVertex++){
			vertCoords[aVertex][0] = mpas1.xVertex[mpas1.vertexIndex[mpas1.verticesOnCell[(iCell)*mpas1.maxEdges+aVertex]]];
			vertCoords[aVertex][1] = mpas1.yVertex[mpas1.vertexIndex[mpas1.verticesOnCell[(iCell)*mpas1.maxEdges+aVertex]]];
			vertCoords[aVertex][2] = mpas1.zVertex[mpas1.vertexIndex[mpas1.verticesOnCell[(iCell)*mpas1.maxEdges+aVertex]]];

			// set areaB[aVertex] here later
			// dprint("iCell %d, vertexIndex[] %d", iCell, mpas1.vertexIndex[mpas1.verticesOnCell[(iCell)*mpas1.maxEdges+aVertex]]);
			// fprintf(stderr, "%f %f %f ", vertCoords[aVertex][0], vertCoords[aVertex][1], vertCoords[aVertex][2]);
		}

		// exit(0);

		// initialize velocities to 0
	    particleVelocity = 0;
	    particleVelocityVert = 0;

	    // general interpolation for the velocity field

	    for (int aTimeLevel=0; aTimeLevel<timeInterpOrder; aTimeLevel++){

	    	// get uVertexVelocity, vVertexVelocity, wVertexVelocity for current time: done implicitly by update_velocity_vectors()

	    	particle_vertical_treatment(nCellVertices, &mpas1.verticesOnCell[iCell*mpas1.maxEdges], uVertexVelocities[aTimeLevel], vVertexVelocities[aTimeLevel], wVertexVelocities[aTimeLevel], mpas1, uvCell, iLevel, zSubStep, &zMid_cur[iCell*mpas1.nVertLevels], &zTop_cur[iCell*mpas1.nVertLevels], &vertVelocityTop_cur[iCell*mpas1.nVertLevelsP1], verticalVelocityInterp);



	    	 particleVelocityVert = particleVelocityVert + timeCoeff[aTimeLevel] * verticalVelocityInterp;


	    	 particleVelocity = particleVelocity + timeCoeff[aTimeLevel]* particle_horizontal_interpolation(mpas1, nCellVertices, vertCoords, xSubStep, uvCell, areaB);

	    }

	    // dprint("Exiting"); exit(0);

	}

	Array3d pathline::particle_horizontal_interpolation(mpas_io &mpas1, const int nCellVertices, const double vertCoords[][3], const Array3d &pointInterp, const double uvCell[][3], double *areaB){




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

	double pathline::wachspress_interpolate(const double *lambda, const double phi[][3], const int idx, const int nCellVertices ){

		double sum = 0;

		for (int i=0; i<nCellVertices; i++)
			sum += lambda[i]*phi[i][idx];

		return sum;
	}


	void pathline::wachspress_coordinates(const mpas_io &mpas1, const int nVertices, const double vertCoords[][3], const Array3d &pointInterp, double *areaB, double *wachspress){


		double areaA[nVertices], wach[nVertices];

		double p[3] = {pointInterp(0), pointInterp(1), pointInterp(2)};
		int im1, i0, ip1;
		for (int i=0; i<nVertices; i++){
			im1 = (nVertices + i - 2)% nVertices ;
			i0  = (nVertices + i - 1)% nVertices ;
			ip1 = (nVertices + i ) % nVertices ;

			double a[3] = {vertCoords[im1][0], vertCoords[im1][1], vertCoords[im1][2]};
			double b[3] = {vertCoords[i0][0], vertCoords[i0][1], vertCoords[i0][2]};
			double c[3] = {vertCoords[ip1][0], vertCoords[ip1][1], vertCoords[ip1][2]};

			areaB[i] = triangle_signed_area(a, b, c, mpas1.radius);
			areaA[i] = triangle_signed_area(p, b, c, mpas1.radius);

		// dprint("areaA[i] %.20g areaB[i] %.20g", areaA[i], areaB[i]);
		// dprint("area %.20g", areaB[i]);
		}

		/* For each vertex compute wachspress coordinate */
		for (int i=0; i<nVertices; i++){
			 wach[i] = areaB[i];
	        for (int j = i+1; j< i + nVertices - 1; j++){

	        	
	            i0  = (nVertices + j)% nVertices;
	            // dprint("%d %d", i, i0);
	           //accumulate products for A_ij subareas
	            wach[i] = wach[i] * areaA[i0];
	        }
	        // dprint(" ");

		}

		double wach_total = 0;
		for (int i=0; i<nVertices; i++){
			wach_total += wach[i];
		}
	// dprint("wach_total %20f", wach_total);

	/* Compute lambda */

		// std::vector<double> wachspress(nVertices);
		wachspress[nVertices-1] = wach[0]/wach_total;
		for (int i=1; i<nVertices;i++){
			wachspress[i-1] = wach[i]/wach_total;
		}

	}




	void pathline::get_bounding_indices(int &iHigh, int &iLow, const double phiInterp, const double *phiVals, const int iLevel, const int nVertLevels){

		if (phiInterp > phiVals[iLevel]){
			iHigh = iLevel + 1;
			iLow = iLevel;

			// todo:001 should not come here!!!

			// dprint("Shouldn't be here");
			// assert(0);

		}else{

			iHigh = iLevel;
			iLow = iLevel + 1;
		}
			
		 // check to make sure points are in range
	     // optimization point: smarter algorithm won't have to call this ever
	     if(!((phiInterp <= phiVals[iHigh]) && (phiInterp >= phiVals[iLow]))){
	     	// todo:002
	     	// assert(0);
	     	get_bounding_indices_brute_force(nVertLevels, phiInterp, phiVals, iLow, iHigh);
	     }


	}

	void pathline::get_bounding_indices_brute_force(const int nVertLevels, const double phiInterp, const double *phiVals, int &iHigh, int &iLow){

		// make no assumptions
		for (int aLevel=0; aLevel<nVertLevels-1; ++aLevel){
			if (phiVals[aLevel]<=phiInterp && phiInterp <= phiVals[aLevel+1]){
				iLow  = aLevel;
         		iHigh = aLevel+1;
         		break;
			}else if (phiVals[aLevel+1]<=phiInterp && phiInterp <=phiVals[aLevel]){
				iLow  = aLevel+1;
         		iHigh = aLevel;
         		break;
			}

		}
		// assert(phiVals[iLow] <= phiInterp);
		// assert(phiInterp <= phiVals[iHigh]);
	}


	void pathline::particle_vertical_treatment(const int nCellVertices, const int *verticesOnCell, const double *uVertexVelocity, const double *vVertexVelocity, const double *wVertexVelocity, mpas_io &mpas1,double uvCell[][3], const int iLevel, const double zLoc, const double *zMid, const double *zTop, const double *zVertVelocityTop, double &verticalVelocityInterp){


		// vertical treatment = passiveFloat

		// interpolate vertical velocity
		verticalVelocityInterp = interp_vert_velocity_to_zlevel(mpas1, iLevel, zLoc, zTop, zVertVelocityTop);



		// interpolate the horizontal velocity based on z-levels
		interp_nodal_vectors(nCellVertices, verticesOnCell, iLevel, mpas1.nVertLevels, zLoc, zMid, uVertexVelocity, vVertexVelocity, wVertexVelocity, uvCell);

		// dprint("uvCell %f %f %f, %f %f %f", uvCell[0][0], uvCell[0][1], uvCell[0][2], uvCell[1][0], uvCell[1][1], uvCell[1][2]);

	}

	void pathline::interp_nodal_vectors(const int nCellVertices, const int *verticesOnCell, const int iLevel, int nVertLevels, const double phiInterp, const double *phiVals, const double *uVertexVelocity, const double *vVertexVelocity, const double *wVertexVelocity, double uvCell[][3]){

		if (iLevel<1){
			if (iLevel == 0){
				for (int aVertex=0; aVertex<nCellVertices; aVertex++){
					uvCell[aVertex][0] = uVertexVelocity[0]; // 0 = maxloc(phiVals(1:nVertLevels),1)
					uvCell[aVertex][1] = vVertexVelocity[0]; // 0 = maxloc(phiVals(1:nVertLevels),1)
					uvCell[aVertex][2] = wVertexVelocity[0]; // 0 = maxloc(phiVals(1:nVertLevels),1)
				}

			}else if(iLevel==-1){
				for (int aVertex=0; aVertex<nCellVertices; aVertex++){
					uvCell[aVertex][0] = uVertexVelocity[nVertLevels-1]; // nVertLevels-1 = minloc(..)
					uvCell[aVertex][1] = vVertexVelocity[nVertLevels-1]; // nVertLevels-1 = minloc(..)
					uvCell[aVertex][2] = wVertexVelocity[nVertLevels-1]; // nVertLevels-1 = minloc(..)
				}
			}
		}

		int iHigh, iLow;
		const double eps = 10e-14;
		double alpha = 0;
		int theVertex;

		
		// dprint("iLevel %d, iHigh %d, iLow %d, phiInterp %f", iLevel, iHigh, iLow, phiInterp);
		get_bounding_indices(iHigh, iLow, phiInterp, phiVals, iLevel, nVertLevels);

		// dprint("iLevel %d, iHigh %d, iLow %d, phiInterp %f, phiVals[iLo iHi] %f %f", iLevel, iHigh, iLow, phiInterp, phiVals[iLow], phiVals[iHigh]);
		// dprint("uVertexVelocity[0,1,2] %f %f %f", uVertexVelocity[theVertex*nVertLevels], uVertexVelocity[theVertex*nVertLevels+1], uVertexVelocity[theVertex*nVertLevels+2]);

		// dprint("iLevel %d, iHigh %d, iLow %d, phiInterp %f, phiVals[iLevel] %d", iLevel, iHigh, iLow, phiInterp, 0);

		// interpolate to vertical level now
		if (abs(phiVals[iHigh] - phiVals[iLow])<eps){
			alpha = -0.5;
		}else{
			 alpha = (phiInterp - phiVals[iLow])/(phiVals[iHigh] - phiVals[iLow]);
		}

		// dprint("val %f", *uVertexVelocity);
		// interpolate to the vertical level
		for (int aVertex=0; aVertex<nCellVertices; ++aVertex){
			theVertex = verticesOnCell[aVertex]-1;
			uvCell[aVertex][0] = alpha * uVertexVelocity[theVertex*nVertLevels + iHigh] + (1 - alpha) * uVertexVelocity[theVertex*nVertLevels + iLow];
			uvCell[aVertex][1] = alpha * vVertexVelocity[theVertex*nVertLevels + iHigh] + (1 - alpha) * vVertexVelocity[theVertex*nVertLevels + iLow];
			uvCell[aVertex][2] = alpha * wVertexVelocity[theVertex*nVertLevels + iHigh] + (1 - alpha) * wVertexVelocity[theVertex*nVertLevels + iLow];

			// dprint("theVertex %d, iHigh %d, uVertexVelocity[theVertex*nVertLevels + iHigh] %f", theVertex, iHigh, uVertexVelocity[theVertex*nVertLevels + iHigh]);

		}

	}

	double pathline::interp_vert_velocity_to_zlevel(mpas_io &mpas1, const int iLevel, const double zSubStep, const double *zTop, const double *vertVelocityTop){


		double alpha=0, interp_vert_velocity_to_zlevel;
		// dprint("iLevel %d, alpha %f, vertVelocityTop[iLevel] %f", iLevel, alpha, vertVelocityTop[0]);

		if (iLevel < 1){
			if (iLevel==0){
				interp_vert_velocity_to_zlevel = vertVelocityTop[0];
			}else if (iLevel==-1){
				interp_vert_velocity_to_zlevel = vertVelocityTop[mpas1.nVertLevelsP1-1];
			}
		}else{

			// interpolate the vertical velocity [assumes zTop(iLevel+1) <= zSubStep <= zTop(iLevel)]
			if (zTop[iLevel+1] <= zSubStep && zSubStep <= zTop[iLevel]){
         		alpha = (zSubStep - zTop[iLevel+1])/(zTop[iLevel]- zTop[iLevel+1]);
       		}else if (zTop[iLevel] <= zSubStep && zSubStep <= zTop[iLevel+1]){
         		alpha = (zSubStep - zTop[iLevel])/(zTop[iLevel+1]- zTop[iLevel]);
     		}

			interp_vert_velocity_to_zlevel = alpha*vertVelocityTop[iLevel] + (1-alpha)*vertVelocityTop[iLevel+1];
		}
		
		// dprint("interp_vert_velocity_to_zlevel %f", interp_vert_velocity_to_zlevel);
		return interp_vert_velocity_to_zlevel;
		
	}



void pathline::update_velocity_vectors(mpas_io &mpas1, int frame_no){

	uVertexVelocities[0] = &mpas1.velocityXv[(frame_no+1)%2][0];
	vVertexVelocities[0] = &mpas1.velocityYv[(frame_no+1)%2][0];
	wVertexVelocities[0] = &mpas1.velocityZv[(frame_no+1)%2][0];

	uVertexVelocities[1] = &mpas1.velocityXv[frame_no%2][0];
	vVertexVelocities[1] = &mpas1.velocityYv[frame_no%2][0];
	wVertexVelocities[1] = &mpas1.velocityZv[frame_no%2][0];

	zMid_cur = &mpas1.zMid[(frame_no+1)%2][0];
	zTop_cur = &mpas1.zTop[(frame_no+1)%2][0];
	vertVelocityTop_cur = &mpas1.vertVelocityTop[(frame_no+1)%2][0];
}



