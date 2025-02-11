#ifndef COMPLEX_ARRAY_H
#define COMPLEX_ARRAY_H

#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace ModuleBase
{

class ComplexArray
{
public:
//	const static int datasize;	// = sizeof(std::complex);
	std::complex<double> *ptr=nullptr; // data array
	
	ComplexArray(const int bnd1=0, const int bnd2=1, const int bnd3=1, const int bnd4=1);

	~ComplexArray();

	void freemem();

	void create(const int bnd1=0, const int bnd2=1, const int bnd3=1, const int bnd4=1);

	ComplexArray(const ComplexArray &cd);
	ComplexArray(ComplexArray &&cd);
	ComplexArray& operator=(ComplexArray &&cd);

	// get and release scratch space
//  void get_temp(int);
//  void release_temp();

	ComplexArray &operator=(const ComplexArray &cd);
	inline void operator=(std::complex <double> c);
//  inline std::complex < double>  &operator()(int i, int j, int k)const
//{return d[(i * bound2 + j) * bound3 +k];}
//  inline void operator=(std::complex < double> c);

	ComplexArray operator+(const ComplexArray &cd);
	void operator+=(const ComplexArray &cd);
	ComplexArray operator-(const ComplexArray &cd);
	void operator-=(const ComplexArray &cd);
	ComplexArray operator*(const double r);
	ComplexArray operator*(const std::complex < double> c);
	void operator*=(const double r);
	void operator*=(const std::complex < double> c);

	void operator*=(const ComplexArray &in);

	// subscript operator
	std::complex < double> &operator()
		(const int ind1=0, const int ind2=0, const int ind3=0, const int ind4=0);
//  std::complex < double> &operator()(int, int, int, int, int);

	const std::complex < double> &operator()
		(const int ind1=0, const int ind2=0, const int ind3=0, const int ind4=0)const;
//  const std::complex < double> &operator()(int, int, int, int, int)const;

	void zero_out(void);
	void negate(void);
	void randomize(void);//uniform distribution

//  void write(char *fname);
//  void write(FILE *fp);
//  void writea(char *fname);
//  void read(char *fname);
	void print();
	
	int getBound1()const{ return bound1; }
	int getBound2()const{ return bound2; }
	int getBound3()const{ return bound3; }
	int getBound4()const{ return bound4; }
	int getSize()const{ return bound1*bound2*bound3*bound4; }

private:
	int bound1, bound2, bound3, bound4;	
	void init(const int size);
};

ComplexArray operator*(double r, const ComplexArray &cd);
ComplexArray operator*(std::complex < double> c, const ComplexArray &cd);

double abs2(const ComplexArray &cd);

void
add_scale_abs2(const std::complex < double> &c, const ComplexArray & in,
               ComplexArray &out);

std::complex < double> dot(const ComplexArray &cd1, const ComplexArray &cd2);
void scale_accumulate(double r, const ComplexArray &cd1, ComplexArray &cd2);
void scale_accumulate(std::complex < double> c, const ComplexArray &cd1, ComplexArray &cd2);
void scaled_sum(double r1, const ComplexArray &cd1,
                double r2, const ComplexArray &cd2,
                ComplexArray &cd3);

void scaled_sum(double r1, ComplexArray &cd1,
                double r2, ComplexArray &cd2,
                ComplexArray &cd3);

void point_mult(ComplexArray &a1, ComplexArray &in2, ComplexArray &out);

// set elements of u as zero which u is 1_d std::complex array
template <class T>
void zeros(std::complex <T> *u, int n)
{
	if (n == 0 || u == 0)
	{
		std::cout << "\n error in zeros(),n or u = 0";
		return;
	}

	for (int i = 0;i < n;i++)
	{
		u[i] = std::complex <T> (0.0, 0.0);
	}
}
}

#endif // COMPLEX_ARRAY_H
