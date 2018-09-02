#ifndef FLOW_H
#define FLOW_H

#include <list>
#include <string>
#include <vector>

class flow
{



public:
    int num_lines;
    std::vector<int> fl_lengths; // stores lengths of flowlines


    std::list<double> flow_lines; // time step, line, dimensions


    flow(int num_lines);
    ~flow();

    void write_flow_to_vtk(std::string &filename);

};

#endif
