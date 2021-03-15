#ifndef BLOCK_H
#define BLOCK_H

#include "mpas_io.h"
#include "ptrace.h"
#include <vector>
#include <string>
#include <set>
#include <diy/mpi.hpp>
// #include <pnetcdf.h>

#include <diy/master.hpp>
#include <Eigen/Dense>


class block: public mpas_io{



public:

	int gid;
	size_t maxEdges = 6;
	// size_t nCells = 7234;

	std::vector<int> gcIdToGid; // global cell index to global block id
	std::vector<int> gVIdxToGid; // global vertex index to global block id
	std::vector<Halo> halo_info; // stores the info on halo cellID and vertexID

	/* To check if inside partition */
	// vector<int> in_partition;

	// vector<int> in_partition_init;

	

	/*particle related */
	std::vector<EndPt> particles;
	int init=0, done=0;
	vector<Segment> segments; // finished segments of particle traces
	std::vector<EndPt> particles_store; // staging area for dequeing particles (incomplete in current frame)
	std::vector<EndPt> particles_continuing; // for particles continuing across epochs/frames in pathlines

	// for debugging
	int framenum;

	std::vector<int> local_gcIds;
	std::vector<int> gcIdToGid_local;
	std::vector<int> gcIdToGid_init; // for original partition
	std::vector<int> local_gcIds_init; // for original partition

	void create_links(const std::string &fname_graph, const std::string &fname_graphpart, std::set<int> &links);

	// void create_links_mpas(const std::string &fname_particles, std::set<int> &links, diy::mpi::communicator &world);
	
	void init_seeds_mpas(std::string &fname_particles, int framenum, int rank);

	void generate_new_particle_file();

	std::vector<Halo> get_halo_info();
	void process_halo_req(Halo &h, int framenum);
	void update_halo_info(Halo &h, int framenum);

	void process_halo_dynamic(int framenum);
	void update_halo_dynamic(Halo &h, int framenum);


	void parallel_write_simstep_segments(diy::mpi::communicator &comm, int framenum);
	void parallel_write_segments(diy::mpi::communicator &comm, int max_steps);

	void create_links_from_gcIdToGid(const std::string &fname_graph, std::set<int> &links);

	void init_seeds_particles(diy::mpi::communicator& world, std::string &fname_particles, int framenum, int seed_rate);
	// void init_partitions();

	void position_data_before_prediction(diy::mpi::communicator &world, int gid, int framenum, std::vector<EndPt>& particles_hold);

	void position_seeds_before_prediction(int pred_percent, std::vector<EndPt> &particles_hold);
	void move_data_to_preadvection_position(std::vector<int> &local_gcIds, std::vector<int> &gcIdToGid_local, std::vector<EndPt>& particles_hold);

	void position_data_and_seeds_after_balance(diy::mpi::communicator &world, std::vector<EndPt> &particles_hold);

	void copy_cell_data_to_endpt(EndPt &cell_cen, int idx);
	void copy_endpt_to_cell_data(EndPt &cell_cen);

	void update_link(const set<int> &links, const diy::Assigner &assigner, diy::Link* link);

	void cell_to_vertex_interpolation();

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