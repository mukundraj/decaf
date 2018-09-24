#include "block_io.h"
#include <cstdio>
#include <vector>
#include <string.h> /* strcpy() */

#include <diy/mpi.hpp>
#include <pnetcdf.h>


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

