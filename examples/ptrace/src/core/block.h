#ifndef BLOCK_H
#define BLOCK_H

#include "mpas_io.h"
#include "ptrace.h"
#include <vector>
#include <string>
#include <set>
#include <diy/mpi.hpp>

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

	void create_links_mpas(const std::string &fname_particles, std::set<int> &links, diy::mpi::communicator &world);
	
	void init_seeds_mpas(std::string &fname_particles);

	void generate_new_particle_file();

	std::vector<Halo> get_halo_info();
	void process_halo_req(Halo &h, int framenum);
	void update_halo_info(Halo &h, int framenum);
	
	void process_halo_dynamic(Halo &h, int framenum);
	void update_halo_dynamic(Halo &h, int framenum);





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