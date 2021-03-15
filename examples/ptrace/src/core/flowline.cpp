#include "flowline.h"
#include "misc.h"
#include "mpas_io.h"
#include "geometry_utils.h"
#include "block.h"


using namespace Eigen;


void flowline::cell_to_vertex_interpolation(block *b, bool pred_frame){

	
	// // swap the t0 and t1, t0 will be overwritten here
	// std::vector<double> tmp;
	// tmp = std::move(b->uVertexVelocities[0]);
	// b->uVertexVelocities[0] = std::move(b->uVertexVelocities[1]);
	// b->uVertexVelocities[1] = std::move(tmp);
	// tmp = std::move(b->vVertexVelocities[0]);
	// b->vVertexVelocities[0] = std::move(b->vVertexVelocities[1]);
	// b->vVertexVelocities[1] = std::move(tmp);
	// tmp = std::move(b->wVertexVelocities[0]);
	// b->wVertexVelocities[0] = std::move(b->wVertexVelocities[1]);
	// b->wVertexVelocities[1] = std::move(tmp);

	std::vector<int> &local_gcids = (pred_frame==true) ? b->local_gcIds_init : b->local_gcIds; 


	b->uVertexVelocities[0] = std::move(b->uVertexVelocities[1]);
	b->vVertexVelocities[0] = std::move(b->vVertexVelocities[1]);
	b->wVertexVelocities[0] = std::move(b->wVertexVelocities[1]);

	b->uVertexVelocities[1].resize(b->uVertexVelocities[0].size());
	b->vVertexVelocities[1].resize(b->vVertexVelocities[0].size());
	b->wVertexVelocities[1].resize(b->wVertexVelocities[0].size());

	std::unordered_set<int> completed_verts;

	// dprint("verticesOnCell size %ld", b->verticesOnCell.size());


	for (size_t i=0; i<local_gcids.size(); i++){
		int cellIdx = local_gcids[i]-1;

		

		// loop over vertices on cellId
		for (int j=0; j<b->maxEdges; j++){
			int vid = b->verticesOnCell[b->maxEdges*cellIdx + j];

			// if (cellIdx == 5470)
			// 		dprint("FOUND 5470, vid %d | %d | size %ld", vid, completed_verts.find(vid)==completed_verts.end(), completed_verts.size());

			if (cellIdx == 2630-1){
				dprint("FOUND 2630, aVertex %d, b->maxEdges %ld", vid, b->maxEdges);
				// flag = 1;
			}
			
			// if vertId in completed_verts, skip remaining in this iteration
			if (completed_verts.find(vid)==completed_verts.end()){
				completed_verts.insert(vid);
				int aVertex = vid - 1;
				int flag=0;

				int flagg = 0;
				if (vid==870 || vid== 869 || vid==12465 || vid==4989 || vid==4994 || vid==12466){
					dprint("vid %d", vid);
					flagg = vid;
				}

				const int vertexDegree = 3;
				
				// get pointVertex: pos of cells surrounding vertex
				double pointVertex[vertexDegree][3];
				for (int aCell=0; aCell < vertexDegree; aCell++){

					int nbrCellIdx = b->cellsOnVertex[vertexDegree*aVertex + aCell] - 1;
					// pointVertex[0][aCell] = b->xCell[b->cellsOnVertex[cellIdx ]];
					// pointVertex[1][aCell] = b->yCell[b->cellsOnVertex[cellIdx ]];
					// pointVertex[2][aCell] = b->zCell[b->cellsOnVertex[cellIdx ]];
					pointVertex[aCell][0] = b->xCell[nbrCellIdx ];
					pointVertex[aCell][1] = b->yCell[nbrCellIdx ];
					pointVertex[aCell][2] = b->zCell[nbrCellIdx ];

					// if (flagg>0)
					// {
					// 	dprint("vid %d, nbrCellIdx %d, (%f %f %f) %f", vid, nbrCellIdx, b->xCell[nbrCellIdx ], b->yCell[nbrCellIdx ], b->zCell[nbrCellIdx ], pointVertex[0][0]);
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
				wachspress_coordinates(b->radius, vertexDegree, pointVertex, pointInterp, areaB, &lambda[0], flag);

				// if (b->gid == 0)
				// if (flag == 1)
				// dprint("lambda %f %f %f, pointInterp %f %f %f | (%f %f %f) (%f %f %f) (%f %f %f) ", lambda[0], lambda[1], lambda[2], pointInterp(0), pointInterp(1), pointInterp(2), pointVertex[0][0], pointVertex[1][0], pointVertex[2][0], pointVertex[0][1], pointVertex[1][1], pointVertex[2][1], pointVertex[0][2], pointVertex[1][2], pointVertex[2][2]);

				// aif (b->gid==0)
				// 	dprint("b->velocitiesX[1][2629] %f", b->velocitiesX[1][2629*b->nVertLevels]);
				
				// interpolate and store
				for (int aLevel=0; aLevel<b->nVertLevels; aLevel++){
					if (b->boundaryVertex[b->nVertLevels*aVertex + aLevel] < 1){
						
						// get cell center velocities
						double ucCell[vertexDegree][3];

						for(int aCell=0; aCell < vertexDegree; aCell++){
							int nbrCellIdx = b->cellsOnVertex[vertexDegree*aVertex + aCell] - 1;
							ucCell[aCell][0] = b->velocitiesX[1][nbrCellIdx * b->nVertLevels + aLevel];
							ucCell[aCell][1] = b->velocitiesY[1][nbrCellIdx * b->nVertLevels + aLevel];
							ucCell[aCell][2] = b->velocitiesZ[1][nbrCellIdx * b->nVertLevels + aLevel];

							
							
						}

						
						// if (flagg > 0 && aLevel == 67){
						// 		dprint("ucCell %f %f %f, gid %d, aLevel %d", ucCell[0][0], ucCell[0][1], ucCell[0][2], b->gid, aLevel);
						// }
						

						// b->uVertexVelocity[b->nVertLevels*aVertex + aLevel] = wachspress_interpolate(lambda, ucCell, 0, vertexDegree);
						// b->vVertexVelocity[b->nVertLevels*aVertex + aLevel] = wachspress_interpolate(lambda, ucCell, 1, vertexDegree);
						// b->wVertexVelocity[b->nVertLevels*aVertex + aLevel] = wachspress_interpolate(lambda, ucCell, 2, vertexDegree);

						b->uVertexVelocities[1][b->nVertLevels*aVertex + aLevel] = wachspress_interpolate(lambda, ucCell, 0, vertexDegree);
						b->vVertexVelocities[1][b->nVertLevels*aVertex + aLevel] = wachspress_interpolate(lambda, ucCell, 1, vertexDegree);
						b->wVertexVelocities[1][b->nVertLevels*aVertex + aLevel] = wachspress_interpolate(lambda, ucCell, 2, vertexDegree);





					}

					if (flagg > 0 && aLevel == 67){
								dprint("aVertex %d, b->uVertexVelocities[1][b->nVertLevels*aVertex + 67] %f %f %f", aVertex, b->uVertexVelocities[1][b->nVertLevels*aVertex + 67], b->vVertexVelocities[1][b->nVertLevels*aVertex + 67], b->wVertexVelocities[1][b->nVertLevels*aVertex + 67]);
					}
				}


			}

		}

			


	}

	
}


void flowline::get_bounding_indices(int &iHigh, int &iLow, const double phiInterp, const double *phiVals, const int iLevel, const int nVertLevels)
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

void flowline::get_bounding_indices_brute_force(const int nVertLevels, const double phiInterp, const double *phiVals, int &iHigh, int &iLow)
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


Array3d flowline::particle_horizontal_interpolation(double radius, const int nCellVertices, const double vertCoords[][3], const Array3d &pointInterp, const double uvCell[][3], double *areaB)
{

	double lambda[nCellVertices];
	wachspress_coordinates(radius, nCellVertices, vertCoords, pointInterp, areaB, &lambda[0]);

	// dprint("lambda (%f %f %f %f %f %f) vertCoords %f %f %f", lambda[0], lambda[1], lambda[2],lambda[3], lambda[4], lambda[5], vertCoords[0][0], vertCoords[0][1], vertCoords[0][2]);

	// update particle velocities via horizontal interpolation
	Array3d particle_horizontal_interpolation;
	particle_horizontal_interpolation(0) = wachspress_interpolate(lambda, uvCell, 0, nCellVertices);
	particle_horizontal_interpolation(1) = wachspress_interpolate(lambda, uvCell, 1, nCellVertices);
	particle_horizontal_interpolation(2) = wachspress_interpolate(lambda, uvCell, 2, nCellVertices);

	// dprint("particle_horizontal_interpolation (%f %f %f)", particle_horizontal_interpolation(0), particle_horizontal_interpolation(1), particle_horizontal_interpolation(2));

	return particle_horizontal_interpolation;
}


void flowline::particle_vertical_treatment(const int nCellVertices, const int *verticesOnCell, const double *uVertexVelocity, const double *vVertexVelocity, const double *wVertexVelocity, block *b, double uvCell[][3], const int iLevel, const double zLoc, const double *zMid, const double *zTop, const double *zVertVelocityTop, double &verticalVelocityInterp)
{

	// vertical treatment = passiveFloat

	// dprint("here1 zMid[iLevel] %f", zMid[iLevel]);
	// interpolate vertical velocity
	verticalVelocityInterp = interp_vert_velocity_to_zlevel(b, iLevel, zLoc, zTop, zVertVelocityTop);

	// if(b->gid==0)
	// 	dprint("verticesOnCell %d %d %d %d %d %d", verticesOnCell[0], verticesOnCell[1], verticesOnCell[2], verticesOnCell[3], verticesOnCell[4], verticesOnCell[5]);
	// dprint("b->uVertexVelocities[1][870*100] %f, %f", b->uVertexVelocities[1][870*100], uVertexVelocity[870*100] );

	

	// interpolate the horizontal velocity based on z-levels
	interp_nodal_vectors(b,  nCellVertices, verticesOnCell, iLevel, b->nVertLevels, zLoc, zMid, uVertexVelocity, vVertexVelocity, wVertexVelocity, uvCell);

	// dprint("here3");

	// dprint("uvCell %f %f %f, %f %f %f", uvCell[0][0], uvCell[0][1], uvCell[0][2], uvCell[1][0], uvCell[1][1], uvCell[1][2]);

	if (verticesOnCell[0] - 1==6464){
		// dprint("vertVelInterp %f iLevel %d vertVelocityTop[iLevel] %f", verticalVelocityInterp, iLevel, zVertVelocityTop[iLevel]);
	}
}

void flowline::interp_nodal_vectors(block *b, const int nCellVertices, const int *verticesOnCell, const int iLevel, int nVertLevels, const double phiInterp, const double *phiVals, const double *uVertexVelocity, const double *vVertexVelocity, const double *wVertexVelocity, double uvCell[][3])
{	
	// if (b->gid==0)
	// 	dprint("uVertexVelocity[0] %f", uVertexVelocity[0]);
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
				uvCell[aVertex][0] = uVertexVelocity[theVertex * nVertLevels + nVertLevels + 0]; // nVertLevels-1 = minloc(..)
				uvCell[aVertex][1] = vVertexVelocity[theVertex * nVertLevels + nVertLevels + 0]; // nVertLevels-1 = minloc(..)
				uvCell[aVertex][2] = wVertexVelocity[theVertex * nVertLevels + nVertLevels + 0]; // nVertLevels-1 = minloc(..)
			}
		}
	}

	int iHigh, iLow;
	const double eps = 1e-14;
	double alpha = 0;

	// if (b->gid==0){
	// 	dprint("uvCell %f %f %f", uvCell[0][0], uvCell[0][0], uvCell[0][0]);
	// }

	if (b->gid==0)
			dprint("b->uVertexVelocities[1][869*100] %f %f", b->uVertexVelocities[1][869*100], uVertexVelocity[869*100]);
	

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
	// if (b->gid == 0)
	// 	dprint("uvCell %f %f %f %f %f %f, IH %d IL %d", uvCell[0][0], uvCell[1][0], uvCell[2][0], uvCell[3][0], uvCell[4][0], uvCell[5][0], iHigh, iLow);

	// theVertex = verticesOnCell[0] - 1;
	// int theVertex2 = verticesOnCell[1] - 1;
	// int theVertex3 = verticesOnCell[2] - 1;
	// if (theVertex==6464)
	// 	dprint("uvCell %f %f %f", uvCell[0][0], uvCell[1][0], uvCell[2][0]);
}

double flowline::interp_vert_velocity_to_zlevel(block *b, const int iLevel, const double zSubStep, const double *zTop, const double *vertVelocityTop)
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

void flowline::get_validated_cell_id(const block &mpas1, Array3d &xSubStep, int &iCell, int &nCellVertices)
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

void flowline::velocity_time_interpolation(const int timeInterpOrder,
										   const double *timeCoeff,
										   const int iCell,
										   const int iLevel,
										   const int nCellVertices,
										   const double zSubStep,
										   block *b,
										   const Eigen::Array3d &xSubStep,
										   Eigen::Array3d &particleVelocity,
										   double &particleVelocityVert)
{

	double vertCoords[nCellVertices][3];
	double uvCell[nCellVertices][3];
	double verticalVelocityInterp;
	double areaB[nCellVertices];

	// if (b->gid==0)
	// 	dprint("iCell %d, iLevel %d", iCell, iLevel);

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

	// if (b->gid==0)
	// 	dprint("iCell %d", iCell);

		

	// general interpolation for the velocity field
	for (int aTimeLevel = 0; aTimeLevel < timeInterpOrder; aTimeLevel++)
	{

		// get uVertexVelocity, vVertexVelocity, wVertexVelocity for current time: done implicitly by update_velocity_vectors()

		// if (iCell==1000)
		// 	dprint("VONCELL %d %d %d", b->verticesOnCell[iCell * b->maxEdges], b->verticesOnCell[iCell * b->maxEdges +1], b->verticesOnCell[iCell * b->maxEdges +2]);

		// dprint("testval %f", mpas1.zMid[iCell * b->nVertLevels + iLevel]);

		int invertedTimeLevel = (aTimeLevel == 0) ? 1 : 0;

		particle_vertical_treatment(nCellVertices, &b->verticesOnCell[iCell * b->maxEdges], &b->uVertexVelocities[invertedTimeLevel][0], &b->vVertexVelocities[invertedTimeLevel][0], &b->wVertexVelocities[invertedTimeLevel][0], b, uvCell, iLevel, zSubStep, &b->zMid[iCell * b->nVertLevels], &b->zTop[iCell * b->nVertLevels], &b->vertVelocityTop[iCell * b->nVertLevelsP1], verticalVelocityInterp);

		// dprint("b->uVertexVelocity size %ld, uVertVels %f %f %f %f %f %f", b->uVertexVelocity.size(), b->uVertexVelocity[b->nVertLevels*4894], 
		// b->uVertexVelocity[b->nVertLevels*14825], 
		// b->uVertexVelocity[b->nVertLevels*7728], 
		// b->uVertexVelocity[b->nVertLevels*7727], 
		// b->uVertexVelocity[b->nVertLevels*13860], 
		// b->uVertexVelocity[b->nVertLevels*4895]);

		particleVelocityVert = particleVelocityVert + timeCoeff[invertedTimeLevel] * verticalVelocityInterp;

		particleVelocity = particleVelocity + timeCoeff[invertedTimeLevel] * particle_horizontal_interpolation(b->radius, nCellVertices, vertCoords, xSubStep, uvCell, areaB);
	}

	// if (b->verticesOnCell[(iCell)*b->maxEdges + 0] - 1==6464){
	// 	dprint("iCell %d | vel %f %f %f | %f", iCell, particleVelocity(0), particleVelocity(1), particleVelocity(2), particleVelocityVert);
	// }


	// dprint("Exiting"); exit(0);
}



double flowline::wachspress_interpolate(const double *lambda, const double phi[][3], const int idx, const int nCellVertices)
{

	double sum = 0;

	for (int i = 0; i < nCellVertices; i++)
		sum += lambda[i] * phi[i][idx];

	return sum;
}

void flowline::wachspress_coordinates(double radius, const int nVertices, const double vertCoords[][3], const Array3d &pointInterp, double *areaB, double *wachspress, int flag)
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

		areaB[i] = triangle_signed_area(a, b, c, radius, flag);
		areaA[i] = triangle_signed_area(p, b, c, radius, flag);
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

bool flowline::in_global_domain(const Pt& p, double cx, double cy, double cz){

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

