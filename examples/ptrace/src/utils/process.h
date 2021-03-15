#ifndef _PROCESS_H
#define _PROCESS_H

#include <diy/master.hpp>
#include "block.h"



void sort_data(diy::Master &master, const diy::Assigner& assigner);

void enqueue_ghost_cells_to_nbrs(block *b, const diy::Master::ProxyWithLink &cp, const diy::Assigner &assigner, int pred_frame);
void deque_ghost_cells_from_nbrs(block *b, const diy::Master::ProxyWithLink &cp);

void ghost_exchange(diy::Master &master, const diy::Assigner& assigner, bool pred_frame);

void send_sampled_seeds_to_origin(diy::Master &master, const diy::Assigner&                assigner);

void interpolate_cell_to_vertex(diy::Master &master, bool pred_frame);

void print_global_num_particles(std::string &str, diy::mpi::communicator& world, diy::Master &master, int framenum);


#endif