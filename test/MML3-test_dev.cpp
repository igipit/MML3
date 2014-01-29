// MML3-example1.cpp : an example on usng MML3
//


#include"MML3-Matrix.h"
#include"MML3-Vector.h"
#include"MML3-StaticArray.h"
#include"MML3-Timer.h"
#include"MML3-Math.h"
#include<iostream>
#include<iomanip>


int main()
{

	// Vector testing
	{
		MML3::Vector<double> a, b(45);
		MML3::Vector<double> c(b);
		MML3::Vector<double> d = { 1, 2, 3, 4, 5, 6, 7 };
		std::cout << "d: " << d;

		d.resize(5);
		for (size_t i = 1; i <= 5; ++i)
			d(i) = 2 * i;
		std::cout << "d: " << d;

		double* p = d.begin();
		double* e = d.end();
		while (p != e)
			std::cout << *p++ << " ";

	}


	
	std::cout << "\n\nfs_array testing" << std::endl;
	{

	MML3::StaticVector<double, 3,1> a;
	a[0] = 0;
	a[1] = 1;
	a[2] = 2;


	std::cout << a << std::endl;

	a(1) = 1;
	a(2) = 2; a(3) = 3;

	a.print(std::cout);
	//std::array<double, 3> cc{ { 2.1, 3.3, 4.3 } };
	MML3::StaticVector<double, 3,1> c({ 2.1, 3.3, 4.3 });
	std::cout << "c: " << c << std::endl;

	c += a;
	std::cout << "c+=a: " << c << std::endl;
	auto d = c + a;
	std::cout << "d=c+a: " << d << std::endl;

	typedef MML3::StaticArray<double, 3, 3> mat;

	mat A, B( 0 );
	mat C={ 1, 2, 3, 4, 5 };
	mat D = { { 1, 2, 3, 4, 5 }, { 6, 7, 8 }, { 9, 10, 11 }, { 1.1, 2.2, 3.3 } };

	std::cout << "A: " << A << std::endl
		<< "B: " << B << std::endl
		<< "C: " << C << std::endl
		<< "D: " << D << std::endl;

	}
	{
		using MML3::StaticArray;
		typedef StaticArray<double, 4, 1> v4x1_t; // 4-vector
		typedef StaticArray<double, 4, 4> v4x4_t; // 4x4-matrix

		StaticArray<double, 4, 6> A(1.0);
		StaticArray<double, 6, 5> B(2.0);
		auto C = product(A, B);

		std::cout << "C: " << C << std::endl;

		auto D = MML3::transpose(C);

		std::cout << "C': " << D << std::endl;




	}

	{

		typedef MML3::StaticVector<double,3>    V3_t; // 3-vector
		typedef MML3::StaticArray<double, 3, 3> M3_t; // 3x3-matrix

		V3_t	x({ 1, 0, 0 }),
			y({ 0, 1, 0 }),
			z({ 0, 0, 1 });

		V3_t a = 0.2*(x-y) + 0.5*(y+z) + 0.6*z;
		std::cout << "a: " << a << std::endl;

		M3_t A(0.0);
		A(1, 1) = A(2, 2) = A(3, 3) = 1;
		std::cout << "A: " << A << std::endl;








	}




return 0;
}

