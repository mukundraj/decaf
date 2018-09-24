#ifndef PTRACE_H
#define PTRACE_H


struct Pt
{
    float coords[3];                         // (x, y, z)
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

    const float& operator [](int i) const { return pt.coords[i]; }
    float& operator [](int i)             { return pt.coords[i]; }

    EndPt()
        {
            pid      = 0;
            sid      = 0;
	    start_time = 0;
        }
    EndPt(struct Segment& s);                // extract the end point of a segment
};

// one segment of a particle trace (trajectory)
struct Segment
{
    int        pid;                          // particle ID
    vector<Pt> pts;                          // points along trace
    int        sid;                          // segment ID of this part of the trace

    Segment()
        {
            pid      = 0;
            sid      = 0;
        }
    Segment(EndPt& p)                        // construct a segment from one point
        {
            pid      = p.pid;
            sid      = p.sid;
            Pt pt    = { p[0], p[1], p[2] };
            pts.push_back(pt);
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
EndPt(Segment& s)                       // extract the end point of a segment
{
    pid = s.pid;
    sid = s.sid;
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
            }
        static
        void load(diy::BinaryBuffer& bb, Segment& x)
            {
                diy::Serialization<int>::load(bb, x.pid);
                diy::Serialization< vector<Pt> >::
                    load(bb, static_cast< vector<Pt>& >(x.pts));
                diy::Serialization<int>::load(bb, x.sid);
            }
    };
}



#endif

