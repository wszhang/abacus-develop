#ifndef COMPLEXMATRIX_H
#define COMPLEXMATRIX_H

#include <complex>

#include "matrix.h"

#ifdef _MCD_CHECK
#include "../src_parallel/mcd.h"
#endif
namespace ModuleBase
{
class ComplexMatrix
{

public:

	int nr=0;
	int nc=0;
	int size=0;
	std::complex<double> *c=nullptr;

	ComplexMatrix(): nr(0), nc(0), size(0), c(nullptr){}
	ComplexMatrix(const int nrows,const int ncols,const bool flag_zero=true);		// Peize Lin add flag_zero 2019-05-13
	ComplexMatrix(const ComplexMatrix &m1);
	ComplexMatrix(ComplexMatrix && m1);						// Peize Lin add 2016-08-05
	explicit ComplexMatrix(const matrix &m);							// Peize Lin add 2017-03-29
	~ComplexMatrix();
	
	void create(const int nrow,const int ncol,const bool flag_zero=true);		// Peize Lin add flag_zero 2019-05-13
	ComplexMatrix& operator=(const ComplexMatrix &m);
	ComplexMatrix& operator=(ComplexMatrix && m);			// Peize Lin add 2016-08-05

	//============
	// Operators
	//============
	std::complex<double> &operator()(const int ir,const int ic)
	{
		assert(ir>=0);	assert(ir<nr);	assert(ic>=0);	assert(ic<nc);
		return c[ir*nc+ic];//mohan modify in-line 2007-10-1
	}
	const std::complex<double> &operator()(const int ir,const int ic)const
	{
		assert(ir>=0);	assert(ir<nr);	assert(ic>=0);	assert(ic<nc);
		return c[ir*nc+ic];//mohan modify in-line 2007-10-13
	}
	
	ComplexMatrix& operator*=(const std::complex<double> &s);
	ComplexMatrix& operator+=(const ComplexMatrix &m);
	ComplexMatrix& operator-=(const ComplexMatrix &m);
	matrix real() const;						// Peize Lin add 2017-03-29
	
	//==================
	// member function:
	//==================
	void zero_out(void);
	void set_as_identity_matrix(void);

	std::ostream & print( std::ostream & os, const double threshold_abs=0.0, const double threshold_imag=0.0 ) const;		// Peize Lin add 2021.09.08
};

ComplexMatrix operator+(const ComplexMatrix &m1,  const ComplexMatrix &m2);
ComplexMatrix operator-(const ComplexMatrix &m1,  const ComplexMatrix &m2);
ComplexMatrix operator*(const ComplexMatrix &m1,  const ComplexMatrix &m2);
ComplexMatrix operator*(const std::complex<double> &s, const ComplexMatrix &m);
ComplexMatrix operator*(const ComplexMatrix &m,   const std::complex<double> &s);
ComplexMatrix operator*(const double &s,          const ComplexMatrix &m);
ComplexMatrix operator*(const ComplexMatrix &m,   const double &s);
	
std::complex<double> trace(const ComplexMatrix &m);

double abs2_row(const ComplexMatrix &m,const int ir);		// mohan add 2008-7-1
double abs2_column(const ComplexMatrix &m,const int ic);	// mohan add 2008-7-1
double abs2(const ComplexMatrix &m);
double abs2(const int nmat, ComplexMatrix **m);

ComplexMatrix transpose(const ComplexMatrix &m, const bool &conjugate);
ComplexMatrix conj(const ComplexMatrix &m);						// Peize Lin add 2019-05-13

void scale_accumulate(
		const std::complex<double> &s, 
		const ComplexMatrix &min, 
		ComplexMatrix &mout);

void scale_accumulate(
		const int &nmat, 
		const std::complex<double> &s, 
		ComplexMatrix **min, 
		ComplexMatrix **mout);

void scaled_sum(
		const std::complex<double> &s1, 
		const ComplexMatrix &m1, 
		const std::complex<double> &s2, 
		const ComplexMatrix &m2, 
		ComplexMatrix &mout);

void scaled_sum(
		const int &nmat,
		const std::complex<double> &s1, 
		ComplexMatrix **m1, 
		const std::complex<double> &s2,
		ComplexMatrix **m2, 
		ComplexMatrix **mout);
}
#endif
