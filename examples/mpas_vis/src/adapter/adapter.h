#ifndef ADAPTER_H
#define ADAPTER_H
#include <mpi.h>
extern "C"
{

	void func1(int n, double x[]);
        void init_decaf();
        void prod(int n, double x[]);
	void producer_terminate();
        void func3();
        void finish_decaf();
        void f_MPI_Comm_c2f(MPI_Fint &comm);
	void decaf_put_int(int n);
        void stage_put_intarray(int i1d[], int n, int id, int frame_num);
        void stage_put_doublearray(double d1d[], int n, int id, int frame_num);
        void decaf_put(int id, int frame_no,int nVertices, int nCells, int nVertLevels, int n_procs);
}

#endif
