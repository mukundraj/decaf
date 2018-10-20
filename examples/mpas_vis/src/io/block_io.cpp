#include "block_io.h"
#include <cstdio>
#include <vector>
#include <string.h> /* strcpy() */

#include <diy/mpi.hpp>
#include <pnetcdf.h>

#include <pblock.h>
#include "misc.h"


block_io::block_io(){

	// this->b = b;

}


block_io::~block_io(){

}

#define ERR {if(err!=NC_NOERR)printf("Error at line=%d: %s\n", __LINE__, ncmpi_strerror(err));}
/*
void block_io::write_cell_centers( int gid, const diy::mpi::communicator &world, std::vector<double> &xCell){

	char filename[128];
	int i, j, rank, nprocs, verbose=1, err;
	int ncid, cmode, varid, dimid;

	MPI_Comm_size(world, &nprocs);
	MPI_Offset start, count;
	MPI_Offset  global_nx = 5*nprocs;

	// std::vector<int> buf{gid, gid, gid, gid, gid};
	int buf[5] = {gid, gid, gid, gid, gid};

	// fprintf(stderr, "inside write cell centers %f\n",xCell[0]);
	strcpy(filename, "testfile.nc");


	fprintf(stderr, "filename before create %s\n", filename);

	//[> create a new file for writing ----------------------------------------<]
	cmode = NC_CLOBBER | NC_64BIT_DATA;
	err = ncmpi_create(world, filename, cmode, MPI_INFO_NULL, &ncid);
	ERR

	// define dimension x
	err = ncmpi_def_dim(ncid, "x", global_nx, &dimid);

	//[> define a 1D variable of integer type <]
	err = ncmpi_def_var(ncid, "xCell", NC_INT, 1, &dimid, &varid);
	ERR

	//[> do not forget to exit define mode <]
	err = ncmpi_enddef(ncid);
	ERR

        //	[> now we are in data mode <]
	start = gid*5;
	count = 5;
	fprintf(stderr, "inside write cell %lld %lld %d %d %d\n", start, count, buf[0], varid, ncid);
	err = ncmpi_put_vara_int_all(ncid, varid, &start, &count, &buf[0]);
	ERR

		err = ncmpi_close(ncid);
	ERR
}
*/

	void block_io::handle_error(int status, int lineno)
{
	    fprintf(stderr, "Error at line %d: %s\n", lineno, ncmpi_strerror(status));
		MPI_Abort(MPI_COMM_WORLD, 1);

}

void block_io::write_particle_traces(int gid, const diy::mpi::communicator &world, PBlock &b, int max_steps){

	int ret, ncfile, nprocs, rank;
	char filename[256], buf[13] = "Hello World\n";
	int data[2];


	

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
	if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

	ret = ncmpi_def_dim(ncfile, "nParticle", b.global_trace_sizes.size(), &dimid_p);
	// ret = ncmpi_def_dim(ncfile, "nParticle", 4, &dimid_p); 
	if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

	pos_dims[0] = dimid_s;
	pos_dims[1] = dimid_p;
	// define variables
	ret = ncmpi_def_var(ncfile, "xPos", NC_DOUBLE, ndims_pos, pos_dims, &varid_xPos);
	if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

	ret = ncmpi_def_var(ncfile, "yPos", NC_DOUBLE, ndims_pos, pos_dims, &varid_yPos);
	if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

	ret = ncmpi_def_var(ncfile, "zPos", NC_DOUBLE, ndims_pos, pos_dims, &varid_zPos);
	if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

	// MPI_Offset start[ndims_pos], count[ndims_pos];

	// end define mode
	ret = ncmpi_enddef(ncfile);
	if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

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
			if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

			for (int j=0;j<b.segments[i].pts.size();j++)
				buffer_test[j][0] = b.segments[i].pts[j].coords[1];
			ret = ncmpi_put_vara_double_all(ncfile, varid_yPos, start, count, &buffer_test[0][0]);
			if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

			for (int j=0;j<b.segments[i].pts.size();j++)
				buffer_test[j][0] = b.segments[i].pts[j].coords[2];
			ret = ncmpi_put_vara_double_all(ncfile, varid_zPos, start, count, &buffer_test[0][0]);
			if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
			
		}else{

			double buffer_test[1][1];

			start[0] = 0; start[1] = 0;
			count[0] = 0; count[1] = 0;
			ret = ncmpi_put_vara_double_all(ncfile, varid_xPos, start, count, &buffer_test[0][0]);
			if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
			ret = ncmpi_put_vara_double_all(ncfile, varid_yPos, start, count, &buffer_test[0][0]);
			if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
			ret = ncmpi_put_vara_double_all(ncfile, varid_zPos, start, count, &buffer_test[0][0]);
			if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
			
		}
		// dprint("written %d %d", i, b.gid);

	}


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
		size_t nVertices = b.xVertex.size();
		size_t nVertexAllLayers = 100 * nVertices;

		// create buffer arrays
		double *xVertex, *xVelocityXv;
		xVertex = new double [nVertices];
		xVelocityXv = new double [nVertexAllLayers];




		ret = ncmpi_open(MPI_COMM_SELF, filename,
			NC_WRITE|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncfile);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		
		ret = ncmpi_redef(ncfile); /* enter define mode */
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		ret = ncmpi_def_dim(ncfile, "nVertices", nVertices, &dimid_nVertices);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		
		ret = ncmpi_def_dim(ncfile, "nVertexAllLayers", nVertexAllLayers, &dimid_nVertexAllLayers);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		const size_t start_vertices[1] = {0}, size_vertices[1] = {nVertices};
		ret = ncmpi_def_var(ncfile, "xVertex", NC_DOUBLE, 1, &dimid_nVertices, &varid_xVertex);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		ret = ncmpi_def_var(ncfile, "yVertex", NC_DOUBLE, 1, &dimid_nVertices, &varid_yVertex);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		ret = ncmpi_def_var(ncfile, "zVertex", NC_DOUBLE, 1, &dimid_nVertices, &varid_zVertex);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);


		
		const size_t start_velocityV[1] = {0}, size_velocityV[1] = {nVertexAllLayers};
		ret = ncmpi_def_var(ncfile, "velocityVx", NC_DOUBLE, 1, &dimid_nVertexAllLayers, &varid_velocityXv);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		ret = ncmpi_def_var(ncfile, "velocityVy", NC_DOUBLE, 1, &dimid_nVertexAllLayers, &varid_velocityYv);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		ret = ncmpi_def_var(ncfile, "velocityVz", NC_DOUBLE, 1, &dimid_nVertexAllLayers, &varid_velocityZv);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		ret = ncmpi_def_var(ncfile, "zTopVertex", NC_DOUBLE, 1, &dimid_nVertexAllLayers, &varid_zTopVertex);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		
		//	Next define trace_sizes global
		ret = ncmpi_def_var(ncfile, "trace_sizes_global", NC_INT, 1, &dimid_p, &varid_trace_sizes_global);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		


		// end define mode
		ret = ncmpi_enddef(ncfile);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		
		

		MPI_Offset start[1], count[1];
		start[0]=0, count[0]=b.global_trace_sizes.size();
		ret = ncmpi_put_vara_int_all(ncfile, varid_trace_sizes_global, start, count, trace_sizes_global);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);	

		// xVertex copy to buffer arrays and write variable 
		for (size_t i=0;i<nVertices;i++){
			xVertex[i] = b.xVertex[i];
		}
		for (size_t i=0;i<nVertexAllLayers;i++){
			xVelocityXv[i] = b.velocityXv[i];
		}

		start[0]=0, count[0]=nVertices;
		ret = ncmpi_put_vara_double_all(ncfile, varid_xVertex, start, count, xVertex);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		start[0]=0, count[0]=nVertexAllLayers;
		ret = ncmpi_put_vara_double_all(ncfile, varid_velocityXv, start, count, xVelocityXv);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		// yVertex copy to buffer arrays and write variable
		for (size_t i=0;i<nVertices;i++){
			xVertex[i] = b.yVertex[i];
		}
		for (size_t i=0;i<nVertexAllLayers;i++){
			xVelocityXv[i] = b.velocityYv[i];
		}

		start[0]=0, count[0]=nVertices;
		ret = ncmpi_put_vara_double_all(ncfile, varid_yVertex, start, count, xVertex);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		start[0]=0, count[0]=nVertexAllLayers;
		ret = ncmpi_put_vara_double_all(ncfile, varid_velocityYv, start, count, xVelocityXv);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		// zVertex copy to buffer arrays and write variable
		for (size_t i=0;i<nVertices;i++){
			xVertex[i] = b.zVertex[i];
		}
		for (size_t i=0;i<nVertexAllLayers;i++){
			xVelocityXv[i] = b.velocityZv[i];
		}

		start[0]=0, count[0]=nVertices;
		ret = ncmpi_put_vara_double_all(ncfile, varid_zVertex, start, count, xVertex);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		start[0]=0, count[0]=nVertexAllLayers;
		ret = ncmpi_put_vara_double_all(ncfile, varid_velocityZv, start, count, xVelocityXv);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);


		for (size_t i=0;i<nVertexAllLayers;i++){
			xVelocityXv[i] = b.zTopVertex[i];
		}
		start[0]=0, count[0]=nVertexAllLayers;
		ret = ncmpi_put_vara_double_all(ncfile, varid_zTopVertex, start, count, xVelocityXv);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		



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
		size_t nVertices = b.xVertex.size();
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
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		
		ret = ncmpi_def_dim(ncfile, "nVertexAllLayers", nVertexAllLayers, &dimid_nVertexAllLayers);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		const size_t start_vertices[1] = {0}, size_vertices[1] = {nVertices};
		ret = ncmpi_def_var(ncfile, "xVertex", NC_DOUBLE, 1, &dimid_nVertices, &varid_xVertex);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		ret = ncmpi_def_var(ncfile, "yVertex", NC_DOUBLE, 1, &dimid_nVertices, &varid_yVertex);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		ret = ncmpi_def_var(ncfile, "zVertex", NC_DOUBLE, 1, &dimid_nVertices, &varid_zVertex);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);


		
		const size_t start_velocityV[1] = {0}, size_velocityV[1] = {nVertexAllLayers};
		ret = ncmpi_def_var(ncfile, "velocityVx", NC_DOUBLE, 1, &dimid_nVertexAllLayers, &varid_velocityXv);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		ret = ncmpi_def_var(ncfile, "velocityVy", NC_DOUBLE, 1, &dimid_nVertexAllLayers, &varid_velocityYv);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		ret = ncmpi_def_var(ncfile, "velocityVz", NC_DOUBLE, 1, &dimid_nVertexAllLayers, &varid_velocityZv);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		ret = ncmpi_def_var(ncfile, "zTopVertex", NC_DOUBLE, 1, &dimid_nVertexAllLayers, &varid_zTopVertex);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		
		// //	Next define trace_sizes global
		// ret = ncmpi_def_var(ncfile, "trace_sizes_global", NC_INT, 1, &dimid_p, &varid_trace_sizes_global);
		// if (ret != NC_NOERR) handle_error(ret, __LINE__);

		


		// end define mode
		ret = ncmpi_enddef(ncfile);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		
		

		MPI_Offset start[1], count[1];
		// start[0]=0, count[0]=b.global_trace_sizes.size();
		// ret = ncmpi_put_vara_int_all(ncfile, varid_trace_sizes_global, start, count, trace_sizes_global);
		// if (ret != NC_NOERR) handle_error(ret, __LINE__);	

		// xVertex copy to buffer arrays and write variable 
		for (size_t i=0;i<nVertices;i++){
			xVertex[i] = b.xVertex[i];
		}
		for (size_t i=0;i<nVertexAllLayers;i++){
			xVelocityXv[i] = b.velocityXv[i];
		}

		start[0]=0, count[0]=nVertices;
		ret = ncmpi_put_vara_double_all(ncfile, varid_xVertex, start, count, xVertex);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		start[0]=0, count[0]=nVertexAllLayers;
		ret = ncmpi_put_vara_double_all(ncfile, varid_velocityXv, start, count, xVelocityXv);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		// yVertex copy to buffer arrays and write variable
		for (size_t i=0;i<nVertices;i++){
			xVertex[i] = b.yVertex[i];
		}
		for (size_t i=0;i<nVertexAllLayers;i++){
			xVelocityXv[i] = b.velocityYv[i];
		}

		start[0]=0, count[0]=nVertices;
		ret = ncmpi_put_vara_double_all(ncfile, varid_yVertex, start, count, xVertex);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		start[0]=0, count[0]=nVertexAllLayers;
		ret = ncmpi_put_vara_double_all(ncfile, varid_velocityYv, start, count, xVelocityXv);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		// zVertex copy to buffer arrays and write variable
		for (size_t i=0;i<nVertices;i++){
			xVertex[i] = b.zVertex[i];
		}
		for (size_t i=0;i<nVertexAllLayers;i++){
			xVelocityXv[i] = b.velocityZv[i];
		}

		start[0]=0, count[0]=nVertices;
		ret = ncmpi_put_vara_double_all(ncfile, varid_zVertex, start, count, xVertex);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		start[0]=0, count[0]=nVertexAllLayers;
		ret = ncmpi_put_vara_double_all(ncfile, varid_velocityZv, start, count, xVelocityXv);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);


		for (size_t i=0;i<nVertexAllLayers;i++){
			xVelocityXv[i] = b.zTopVertex[i];
		}
		start[0]=0, count[0]=nVertexAllLayers;
		ret = ncmpi_put_vara_double_all(ncfile, varid_zTopVertex, start, count, xVelocityXv);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);
		



		// close file
		ret = ncmpi_close(ncfile);
		if (ret != NC_NOERR) block_io::handle_error(ret, __LINE__);

		// free buffer arrays
		delete [] xVertex;
		delete [] xVelocityXv;


	}
}

// test parallel writing function
	void block_io::write_cell_centers( int gid, const diy::mpi::communicator &world, std::vector<double> &xCell){

		int ret, ncfile, nprocs, rank, dimid, varid1, varid2, ndims=1;
		MPI_Offset start, count=2;
		char filename[256], buf[13] = "Hello World\n";
		int data[2];


		MPI_Comm_rank(world, &rank);
		MPI_Comm_size(world, &nprocs);


		strcpy(filename, "testfile.nc");

		ret = ncmpi_create(world, filename,
				NC_CLOBBER|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncfile);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		ret = ncmpi_def_dim(ncfile, "d1", nprocs*2, &dimid);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		ret = ncmpi_def_var(ncfile, "v1", NC_INT, ndims, &dimid, &varid1);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		ret = ncmpi_def_var(ncfile, "v2", NC_INT, ndims, &dimid, &varid2);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		ret = ncmpi_put_att_text(ncfile, NC_GLOBAL, "string", 13, buf);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		// all processors defined the dimensions, attributes, and variables,
		//      * but here in ncmpi_enddef is the one place where metadata I/O
		//           * happens.  Behind the scenes, rank 0 takes the information and writes
		//                * the netcdf header.  All processes communicate to ensure they have
		//                     * the same (cached) view of the dataset
                //
		ret = ncmpi_enddef(ncfile);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		start=rank*2, count=2, data[0]=rank, data[1]=rank;

		ret = ncmpi_put_vara_int_all(ncfile, varid1, &start, &count, &data[0]);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

		//ret = ncmpi_put_vara_int_all(ncfile, varid2, &start, &count, &data[0]);
		//if (ret != NC_NOERR) handle_error(ret, __LINE__);

		ret = ncmpi_close(ncfile);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

	}

