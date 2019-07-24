#include "ptrace.h"



const double& EndPt::operator [](int i) const { return pt.coords[i]; }
double& EndPt::operator [](int i)             { return pt.coords[i]; }

EndPt::EndPt()
{
	pid      = 0;
	sid      = 0;
	nsteps   = 0;
}

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


 Segment::Segment()
        {
            pid      = 0;
            sid      = 0;
        }
 Segment::Segment(EndPt& p)                        // construct a segment from one point
        {
            pid      = p.pid;
            sid      = p.sid;
            Pt pt    { { p[0], p[1], p[2] } };
            pts.push_back(pt);
        }