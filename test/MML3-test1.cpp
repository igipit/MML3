// MML3-example1.cpp : an example on usng MML3
//


#include"MML3-Matrix.h"
#include"MML3-Vector.h"
#include"MML3-StaticArray.h"
#include"MML3-Timer.h"
#include"MML3-Math.h"
#include<iostream>
#include<iomanip>
#include"MML3-MatrixAlgorithms.h"

double poly(double x)
{
	double y=0.0;

	const int sz=21;
	double a[sz];
	for( int i=0; i!=sz; ++i)
	{
		a[i]= double(i+1);
		y+=a[i]*pow(x,i);
	}
	return y;
}


using MML3::iSet;
typedef MML3::index_type index_t;
typedef MML3::Mat<>::type		Matrix;
typedef MML3::Mat<>::sym_type symMatrix;



int main()
{

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



	try{

		// a simple stopwatch based on std::chrono
		MML3::Timer timer;

		int test_counter(0), test_failure(0);
		std::cout
			<< "******************************************************\n"
			<< "MML3-test1\n"
			<< "******************************************************\n";


		std::cout << "Test " << ++test_counter << " : direct access to the Matrix components\n";
		{
			Matrix K = Matrix::get_random_matrix(6, 7, 1.0, 0.2);

			for (auto i : K.row_range())
			for (auto j : K.col_range(i))
				K(i, j) = i*j;

			double tmp(0);
			for (auto i : K.row_range())
			for (auto j : K.col_range(i))
				tmp += fabs(K(i, j) - i*j);

			if (fabs(tmp) > 1e-15)
			{
				std::cout << "Error\n";
				test_failure++;
			}
			else
				std::cout << "OK\n";
		}

		std::cout << "Test " << ++test_counter << " : SubMatrix access 1\n";
		{
			Matrix K(5, 6);
			Matrix::iset_t a2({ 1, 2, 3 });
			a2.add(4);
			auto k2 = K(a2, a2);
			k2.fill(1 / 3.0);
			double tmp(0);
			for (int i = 1; i <= 4; ++i)
			for (int j = 1; j <= 4; ++j)
				tmp += K(i, j) - 1 / 3.0;

			if (tmp != 0.)
			{
				std::cout << "Errore\n";
				test_failure++;
			}

			else
				std::cout << "OK\n";
		}



		std::cout << "Test " << ++test_counter << " :  SubMatrix access 2\n";
		{
			Matrix			K(5, 6);
			Matrix::iset_t	RS(1, 3);
			Matrix::iset_t	CS(1, 2);

			auto k = K(RS, CS);
			double tmp(0);
			for (size_t i = 1; i <= RS.size(); ++i)
			for (size_t j = 1; j <= CS.size(); ++j)
				tmp += fabs(K(RS(i), CS(j)) - k(i, j));

			if (tmp != 0.)
			{
				std::cout << "Errore\n";
				test_failure++;
			}

			else
				std::cout << "OK\n";
		}


		std::cout << "\nTest " << ++test_counter << " : GEMM matrix product C= A * B\n";
		{
			int sz = 500, sz1 = 478;


			Matrix  C(sz, sz), C1(sz, sz);
			Matrix	A = Matrix::get_random_matrix(sz, sz1, 0.0, 10.0);
			Matrix  B = Matrix::get_random_matrix(sz1, sz, 1.0, 2.0);

			double start = timer.seconds();
			MML3::Algorithm::product(A, B, C);
			double t1 = timer.seconds();
			// Home Made product
			for (auto i = 0; i != sz; ++i)
			{
				for (auto j = 0; j != sz; ++j)
				{
					double tmp = 0;
					for (auto k = 0; k != sz1; ++k)
						tmp += A[i][k] * B[k][j];
					C1[i][j] = tmp;
				}
			}
			double t2 = timer.seconds();
			C -= C1;
			double err = MML3::Algorithm::norm2(C);

			if (err > 1E-15)
			{
				std::cout << "Errore\n";
				test_failure++;
			}
			else{
				std::cout << "OK\n"
					<< "product() time: " << (t1 - start) << std::endl
					<< "hm_prod() time: " << (t2 - t1) << std::endl;
			}
		}



		std::cout << "\nTest " << ++test_counter << " : GEMM matrix product C= A' * B\n";
		{
			int sz = 100, sz1 = 78;


			Matrix  C(sz, sz), C1(sz, sz);
			Matrix	A = Matrix::get_random_matrix(sz1, sz, 0.0, 10.0);
			Matrix  B = Matrix::get_random_matrix(sz1, sz, 1.0, 2.0);

			double start = timer.seconds();
			MML3::Algorithm::product(A, B, C, true);
			double t1 = timer.seconds();
			for (auto i = 0; i != sz; ++i)
			{
				for (auto j = 0; j != sz; ++j)
				{
					double tmp = 0;
					for (auto k = 0; k != sz1; ++k)
						tmp += A[k][i] * B[k][j];
					C1[i][j] = tmp;
				}
			}
			double t2 = timer.seconds();
			C -= C1;
			double err = MML3::Algorithm::norm2(C);
			if (err > 1E-15)
			{
				std::cout << "Errore\n";
				test_failure++;
			}
			else{
				std::cout << "OK\n"
					<< "product() time: " << (t1 - start) << std::endl
					<< "hm_prod() time: " << (t2 - t1) << std::endl;
			}
		}


		std::cout << "\nTest " << ++test_counter << " : GEMM matrix product C= A * B'\n";
		{
			int sz = 100, sz1 = 78;
			Matrix  C(sz, sz), C1(sz, sz);
			Matrix	A = Matrix::get_random_matrix(sz, sz1, 0.0, 10.0);
			Matrix  B = Matrix::get_random_matrix(sz, sz1, 1.0, 2.0);

			/*int sz = 3, sz1 = 2;
			Matrix  C(sz, sz), C1(sz, sz);
			Matrix	A = { { 1, 2 }, { 10, 20 }, { 100, 200 } };
			Matrix  B = { { 1, 1 }, { 2, 3 },{3,3} };*/

			double start = timer.seconds();
			MML3::Algorithm::product(A, B, C, false, true);
			double t1 = timer.seconds();
			for (auto i = 0; i != sz; ++i)
			{
				for (auto j = 0; j != sz; ++j)
				{
					double tmp = 0;
					for (auto k = 0; k != sz1; ++k)
						tmp += A[i][k] * B[j][k];
					C1[i][j] = tmp;
				}
			}
			double t2 = timer.seconds();
			C -= C1;
			double err = MML3::Algorithm::norm2(C) / MML3::Algorithm::norm2(C1);
			if (err > 1E-14)
			{
				std::cout << "Errore: " << err << std::endl;
				test_failure++;
			}
			else{
				std::cout << "OK  (rdiff=" << err<<")\n"
					<< "product() time: " << (t1 - start) << std::endl
					<< "hm_prod() time: " << (t2 - t1) << std::endl;
			}
		}


		std::cout << "\nTest " << ++test_counter << " : GEMM matrix product C= A' * B'\n";
		{
			int sz = 100, sz1 = 78;


			Matrix  C(sz, sz), C1(sz, sz);
			Matrix	A = Matrix::get_random_matrix(sz1, sz, 0.0, 10.0);
			Matrix  B = Matrix::get_random_matrix(sz, sz1, 1.0, 2.0);

			double start = timer.seconds();
			MML3::Algorithm::product(A, B, C, true, true);
			double t1 = timer.seconds();
			for (auto i = 0; i != sz; ++i)
			{
				for (auto j = 0; j != sz; ++j)
				{
					double tmp = 0;
					for (auto k = 0; k != sz1; ++k)
						tmp += A[k][i] * B[j][k];
					C1[i][j] = tmp;
				}
			}
			double t2 = timer.seconds();
			C -= C1;
			double err = MML3::Algorithm::norm2(C);
			if (err > 1E-15)
			{
				std::cout << "Errore\n";
				test_failure++;
			}
			else{
				std::cout << "OK\n"
					<< "product() time: " << (t1 - start) << std::endl
					<< "hm_prod() time: " << (t2 - t1) << std::endl;
			}
		}



		{
			std::cout << "\nTest " << ++test_counter << " : matrix fread,fwrite\n";
			int sz = 20;
			Matrix A, C = Matrix::get_random_matrix(sz - 1, sz, 0.0, 10.0);
			C.fwrite("pippo");
			A.fread("pippo");
			A -= C;
			double err = MML3::Algorithm::norm2(A);
			if (err != 0.0)
			{
				std::cout << "Errore\n";
				test_failure++;
			}
			else
				std::cout << "OK\n";
		}



		{
			std::cout << "\nTest " << ++test_counter << " : matrix imax,imin\n";
			// test imax imin

			std::cout << "Test imax :";
			Matrix A = {  100., -10, 1, 2, 3, 4, 5, 7, 9, 1000, -34 };
			index_t imax1 = MML3::Algorithm::imax(A);

			if (imax1 != 10)
			{
				std::cout << "Errore\n";
				test_failure++;
			}
			else
				std::cout << "OK\n";



			std::cout << "Test imin :";
			A = { 1, 2, 5, 7, -3, -10, -100, -7, -9 };
			index_t imin1 = MML3::Algorithm::imin(A);

			//std::cout << A << std::endl;

			if (imin1 != 7)
			{
				{
					std::cout << "Errore\n";
					test_failure++;
				}
			}
			else
				std::cout << "OK\n";
		}


		{

			std::cout << "******    component access  OP " << std::endl;
			size_t NR = 100;

			const size_t csz = 500;
			size_t csz1 = csz + 1;
			
			Matrix    xcm(csz, csz);
			double    cm1d[csz*csz];

			double acc_i = 0;
			double start0 = timer.now();
			for (size_t rep = 0; rep < NR; ++rep)
			{
				for (size_t i = 1; i != csz1; ++i)
				{
					double * wi = &cm1d[(i - 1)*csz];
					for (size_t j = 1; j != csz1; ++j)
						wi[  j - 1] = j + i;
					acc_i += wi[ i - 1];
				}
			}
			double elt0 = timer.now() - start0;
			//std::cout << " c access [][] time: " << elt0 << std::endl;



			for (size_t rep = 1; rep <= NR; ++rep)
			{
				for (size_t i = 1; i != csz1; ++i)
				{
					for (size_t j = 1; j != csz1; ++j)
						xcm(i, j) = j + i;
					acc_i += xcm(i, i);
				}
			}

			double elt1 = timer.now() - elt0;
			//std::cout << " MML access () time: " << elt1 << std::endl;
			for (size_t rep = 1; rep <= NR; ++rep)
			{
				for (size_t i = 0; i != csz; ++i)
				{
					double * xi = xcm[i];
					for (size_t j = 0; j != csz; ++j)
						xi[j] = j + i;
					acc_i += xcm[i][i];
				}
			}
			double elt2 = timer.now() - elt1 - elt0;
			// testing ranges access
			for (auto i : xcm.row_range())
			for (auto j : xcm.col_range(i))
				acc_i += xcm(i, j);
			double elt3 = timer.now() - elt2;




			//std::cout << " MML access[][]time: " << elt2 << std::endl;

			double nex = NR*double(csz)*csz;
			std::cout
				<< " C access [][] time: " << elt0 / nex << std::endl
				<< " MML access () time: " << elt1 / nex << "  overload: " << (elt1 / elt0) << std::endl
				<< " MML access[][]time: " << elt2 / nex << "  overload: " << (elt2 / elt0) << std::endl
				<< " MML range     time: " << elt3 / nex << "  overload: " << (elt3 / elt0) << std::endl
				<< " ignore: " << acc_i << std::endl;

			for (size_t i = 1; i <= csz; ++i)
			for (size_t j = 1; j <= csz; ++j)
				xcm(i, j) = cm1d[(i - 1)*csz + j - 1];
			std::cout << "ignore : " << xcm(1, 1) << std::endl;
		}






		{
			typedef MML3::Mat<>::type		Matrix;

			Matrix  A = { { 1, 2, 3, 4 }, { 21, 22, 23, 24 }, { 31, 32, 33, 34 } };
			std::cout << "A" << A << std::endl;

			auto At = MML3::transpose(A);
			std::cout << "A'" << At << std::endl;
		}






	}
	catch (std::exception& e)
	{

		std::cout << "Exception: " << e.what() << std::endl;
		return -1;
	}



return 0;
}

