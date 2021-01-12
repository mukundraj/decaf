#ifndef PARTICLES_IO
#define PARTICLES_IO

#include <vector>
#include <string>


struct seed_data{

    std::vector<double> xPos, yPos, zPos; // seed positions
    std::vector<int> seed_cell_id; // to store the cell id of the cells

    size_t nCells, nVertices;

    std::vector<double> xVertex, yVertex, zVertex;
    std::vector<int> verticesOnCell;
    std::vector<int> nEdgesOnCell;

    
};

void load_MPAS_cells_coords(const std::string& filename, seed_data &data);

void load_particle_data(const std::string& filename, seed_data &data);



#endif