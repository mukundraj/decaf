#ifndef BLOCK_IO
#define BLOCK_IO

// #include "pblock.hpp"
#include <vector>
#include <pnetcdf.h>
#include "mpaso.h"

#include <diy/mpi.hpp>

struct PBlock;
class block_io{
	
private:

	PBlock *b;
public:


	block_io();
	~block_io();
	static void handle_error(int status, int lineno);

	void write_cell_centers(int gid, const diy::mpi::communicator &world, std::vector<double> &xCell);
	void write_particle_traces(int gid, const diy::mpi::communicator &world, PBlock &b, int max_steps);


};


#endif
