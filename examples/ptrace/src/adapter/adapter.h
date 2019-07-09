#ifndef ADAPTER_H
#define ADAPTER_H
#include <mpi.h>
extern "C"
{


	void f_MPI_Comm_c2f(MPI_Fint &comm);
	void decaf_put_int_array(int id, int frame_no, int size, int i1d[]);

	void decaf_put_double_array(int id, int frame_no, int size, double d1d[]);
	void init_decaf();
	void finish_decaf();
}

#endif



