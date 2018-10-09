#ifndef MISC_H
#define MISC_H

#include <string>
#include <vector>
#include <set>


std::string split_filename(std::string str);

#define dprint(format, ...) \
    fprintf(stderr, format " %s %d \n", ## __VA_ARGS__, split_filename(__FILE__).c_str(), __LINE__)

    
std::string itos(int i);


// struct ghost_info
// {
// 	int requester, handler;

// 	std::vector<double> cx,cy,cz,			// cell center pos
// 						vel_cx, vel_cy, vel_cz, // cell center velocities
// 						vert_x, vert_y, vert;	// verterx pos
// 	std::vector<int> cell_gid, vert_gid;					// vertex global ids
	
// };


struct ghost_req
{
	int requester, handler;
	std::set<int> ghost_cell_ids;

	~ghost_req();

};

// TODO: remove the constant 100 and replace by nVertLevels
// template<int nVertLevels>
struct ghost_point
{	
	int cgid; // global ids
	double cx, cy, cz; // cell center pos
	double vel_cx[100], vel_cy[100], vel_cz[100], zTop[100]; // cell center velocities and top


	int vert_gids[6];
	double vert_x[6], vert_y[6], vert_z[6]; // verterx posx, posy, posz
	int cellsOnVertex[18]; // neighbor cell ids for the vertices


};



#endif
