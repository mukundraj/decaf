#ifndef PBLOCK_H
#define PBLOCK_H

#include <vector>
// #include <diy/serialization.hpp>

using namespace std;

struct Pt
{
    double coords[4];                         // (x, y, z)
};

// whether a point is inside given bounds


// bool inside(const Pt& pt, const Bounds bounds)
// {
//     for (int i = 0; i < 3; i++)
//         if (pt.coords[i] < bounds.min[i] || pt.coords[i] >= bounds.max[i] - 1)
//             return false;
//     return true;
// }

// one end point of a particle trace segment
struct EndPt
{
    int  pid;                                // particle ID
    Pt   pt;                                 // end pointof the trace
    int  sid;                                // segment ID of this part of the trace
	int step; 				// number of steps from start for this point in trace
	int gpid;    // global pid based on the block of origin 

    const double& operator [](int i) const { return pt.coords[i]; }
    double& operator [](int i)             { return pt.coords[i]; }

    EndPt()
        {
            pid      = 0;
            sid      = 0;
	    step     = 0; // initially all particles are at step zero
	    gpid     = 0; // need to set this based on the block of initialization
        }
    EndPt(struct Segment& s);                // extract the end point of a segment
};

// one segment of a particle trace (trajectory)
struct Segment
{
    int        pid;                          // particle ID
    vector<Pt> pts;                          // points along trace
    int        sid;                          // segment ID of this part of the trace
    int 	start_step; 			// step of first point in this segment with regard to whole trace
    int gpid; // global pid based on the block of origin of the original segment 
    Segment()
        {
            pid      = 0;
            sid      = 0;
        }
    Segment(EndPt& p  )       // construct a segment from one point.
        {
            pid      = p.pid;
            sid      = p.sid;
            Pt pt    = { p[0], p[1], p[2] };
	    start_step = p.step;
            pts.push_back(pt);
	    gpid = p.gpid;
        }

    // whether end point is inside given bounds
    bool inside(const int lb[3], const int ub[3]) const
        {
            for (int i = 0; i < 3; i++)
                if (pts.back().coords[i] < lb[i] || pts.back().coords[i] >= ub[i] - 1)
                    return false;
            return true;
        }
};

// following constructor defined out of line because references Segment, which needed
// to be defined first
EndPt::
EndPt(Segment& s)                       // extract the end point of a segment.  CHECK also pass start_step here?
{
    pid = s.pid;
    sid = s.sid;
    step = s.start_step + s.pts.size();
    gpid = s.gpid;
    pt.coords[0] = s.pts.back().coords[0];
    pt.coords[1] = s.pts.back().coords[1];
    pt.coords[2] = s.pts.back().coords[2];
}

// specialize the serialization of a segment
namespace diy
{
    template<>
    struct Serialization<Segment>
    {
        static
        void save(diy::BinaryBuffer& bb, const Segment& x)
            {
                diy::Serialization<int>::save(bb, x.pid);
                diy::Serialization< vector <Pt> >::
                    save(bb, static_cast< const vector<Pt>& >(x.pts));
                diy::Serialization<int>::save(bb, x.sid);
                diy::Serialization<int>::save(bb, x.start_step);
                diy::Serialization<int>::save(bb, x.gpid);
            }
        static
        void load(diy::BinaryBuffer& bb, Segment& x)
            {
                diy::Serialization<int>::load(bb, x.pid);
                diy::Serialization< vector<Pt> >::
                    load(bb, static_cast< vector<Pt>& >(x.pts));
                diy::Serialization<int>::load(bb, x.sid);
                diy::Serialization<int>::load(bb, x.start_step);
                diy::Serialization<int>::load(bb, x.gpid);
            }
    };
}

struct PBlock
{
    PBlock(): count(0)                   {}

    	int   count;
	int gid;
    	vector<int> values{9,9,9};
	vector<int> indexToCellID; 	//0
	vector<double> xCell; 		//1
	vector<double> yCell; 		//2
	vector<double> zCell; 		//3
	vector<double> velocityX; 	//4
	vector<double> velocityY; 	//5
	vector<double> velocityZ; 	//6
	vector<int> indexToVertexID; 	//7  
	vector<double> xVertex; 	//8
	vector<double> yVertex; 	//9
	vector<double> zVertex; 	//10
	vector<double> zTop; 		//11
	vector<int> cellsOnVertex; 	//12
    vector<int> verticesOnCell; //13

	 vector<Segment> segments; // finished segments of particle traces
	 vector<int> global_trace_sizes; // latest known size (length) of each block particle to id early finishes 
	 vector<int> global_nP; // global view of particles initialized in each block
	 vector<int> global_start_ids; // global start ids used for computing a global particle id
};

namespace diy
{
    template<>
        struct Serialization<PBlock>
	{
        static void save(BinaryBuffer& bb, const PBlock& b)
	{
            diy::save(bb, b.count);
            diy::save(bb, b.values);
        
	}

        static void load(BinaryBuffer& bb, PBlock& b)
	{
            diy::load(bb, b.count);
            diy::load(bb, b.values);
        
	}
    
	};

}

void* create_block()                      { return new PBlock;  }
void  destroy_block(void* b)              { delete static_cast<PBlock*>(b);  }
void  save_block(const void* b,
                 diy::BinaryBuffer& bb)   { 
						diy::save(bb, *static_cast<const PBlock*>(b));
											
						  
					}
void  load_block(void* b,
                 diy::BinaryBuffer& bb)   { diy::load(bb, *static_cast<PBlock*>(b));  }







#endif
