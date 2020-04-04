#ifndef BLOCK_H
#define BLOCK_H

#include "mpas_io.h"
#include "ptrace.h"
#include <vector>
#include <string>
#include <set>
#include <diy/mpi.hpp>
// #include <pnetcdf.h>


class block: public mpas_io{



public:

	int gid;
	size_t maxEdges = 6;
	// size_t nCells = 7234;

	std::vector<int> gcIdxToGid; // global cell index to global block id
	std::vector<int> gVIdxToGid; // global vertex index to global block id
	std::vector<Halo> halo_info; // stores the info on halo cellID and vertexID

	/*particle related */
	std::vector<EndPt> particles;
	int init=0, done=0;
	vector<Segment> segments; // finished segments of particle traces
	std::vector<EndPt> particles_store; // staging area for dequeing particles


	/* To check if inside partition */
	vector<int> in_partition;
	
	void create_links(const std::string &fname_graph, const std::string &fname_graphpart, std::set<int> &links);
	void create_links_mpas(const std::string &fname_particles, std::set<int> &links, diy::mpi::communicator &world);

	void init_seeds_particles(diy::mpi::communicator& world, std::string &fname_particles, int framenum);
	void init_seeds_mpas(std::string &fname_particles, int framenum, int rank);

	void generate_new_particle_file();

	std::vector<Halo> get_halo_info();
	void process_halo_req(Halo &h, int framenum);
	void update_halo_info(Halo &h, int framenum);

	void process_halo_dynamic(int framenum);
	void update_halo_dynamic(Halo &h, int framenum);


	void parallel_write_simstep_segments(diy::mpi::communicator &comm, int framenum);
	void parallel_write_segments(diy::mpi::communicator &comm, int max_steps);

	void init_partitions();

	

};

namespace diy
{
    template<>
        struct Serialization<block>
    {
        static void save(BinaryBuffer& bb, const block& b)
        {
            diy::save(bb, b.gid);
        }

        static void load(BinaryBuffer& bb, block& b)
        {
            diy::load(bb, b.gid);
        }
    };
}


#endif