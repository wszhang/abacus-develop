#ifndef PZT2S_H
#define PZT2S_H

#include "../src_pw/tools.h"
#include "pdiag_common.h"

void pzt2s(MPI_Comm comm2D,int N_A,int NB,std::complex <double> *A,std::complex <double> *X,
	     LocalMatrix loc_A,int loc_size,std::complex <double> *norm,char uplo);

#endif
