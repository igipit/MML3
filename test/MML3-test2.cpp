// MML3-example1.cpp : an example on usng MML3
//


#include"MML3-Matrix.h"
#include"MML3-Timer.h"
#include"MML3-Math.h"
#include<iostream>
#include<iomanip>
#if defined( MML3_LAPACK ) &&  MML3_LAPACK >0
#include"MML3-MatrixAlgorithms.h"
#endif

#if defined( MML3_LAPACK ) &&  MML3_LAPACK >0
#include"MML3-LinearSolver.h"
#endif

typedef MML3::index_type index_t;

using MML3::Solver::LUSolver;
using MML3::Solver::LLTSolver;
using MML3::Mat;

int main()
{

	// a simple stopwatch based on std::chrono
	MML3::Timer timer;

	int test_counter(0), test_failure(0);


	const size_t sz = 1000, nrhs = 3;
	std::cout << "*********************************************************"
		<< "\nTest on  LU and LLT  factorization and solution  \n"
		<< "size: " << sz << "    rhs: " << nrhs
		<< "\n*********************************************************\n";




	try{

		
		{
			typedef Mat<double>::type		Matrix;
			std::cout << "\nTest " << std::setw(3) << ++test_counter << ":  LU factorization on " << Matrix::name() <<  std::endl;
			
			

			auto A = Matrix::get_random_matrix(sz, sz, 1.0, 2.0);
			for (auto i = 1; i != sz + 1; ++i)
				A(i, i) += 100+i;
			auto B = Matrix::get_random_matrix(sz, nrhs, 1.0, 0.5);
			Matrix X(B), Acp(A);
			timer.start();

			LUSolver<Matrix>::type LU(Acp);

			//MML3::Solver::LU<Matrix> LU(Acp);
			double t0 = timer.seconds();
			if (LU)
				LU.solve(X);
			double t1 = timer.seconds();
			if (!LU)
			{
				std::cout << "Errore\n";
				test_failure++;
			}
			else
			{

				// verify the correctness of the solution 
				Matrix tmp;
				MML3::Algorithm::product(A, X, tmp);
				tmp -= B;
				double err = MML3::Algorithm::norm2(tmp) / MML3::Algorithm::norm2(B);
				if (err > 1e-14)
				{
					std::cout << "Error:" << err << "n";
					test_failure++;
				}
				else
					std::cout << "OK"
					<< "\n factorization time : " << t0
					<< "\n solution time      : " << (t1 - t0) << std::endl;
			}
		}


		{
			
			
			typedef Mat<double>::General::Rectangular::ColMajor::type Matrix;
			typedef Mat<double>::General::Rectangular::ColMajor::type Matrix;
			std::cout << "\nTest " << std::setw(3) << ++test_counter << ":  LU factorization on " << Matrix::name() << std::endl;

			auto A = Matrix::get_random_matrix(sz, sz, 1.0, 2.0);
			for (auto i = 1; i != sz + 1; ++i)
				A(i, i) += 100 + i;
			auto B = Matrix::get_random_matrix(sz, nrhs, 1.0, 0.5);
			Matrix X(B), Acp(A);
			timer.start();
			LUSolver<Matrix>::type  LU(Acp);
			double t0 = timer.seconds();
			if (LU)
				LU.solve(X);
			double t1 = timer.seconds();
			if (!LU)
			{
				std::cout << "Errore\n";
				test_failure++;
			}
			else
			{

				// verify the correctness of the solution 
				Matrix tmp;
				MML3::Algorithm::product(A, X, tmp);
				tmp -= B;
				double err = MML3::Algorithm::norm2(tmp) / MML3::Algorithm::norm2(B);
				if (err > 1e-14)
				{
					std::cout << "Error:" << err << "n";
					test_failure++;
				}
				else
					std::cout << "OK"
								<< "\n factorization time : " << t0
								<< "\n solution time      : " << (t1 - t0) << std::endl;
			}
		}

		
		{
			
			typedef Mat<double>::type Matrix;
			typedef Mat<double>::Symmetric::type symMatrix;

			std::cout << "\nTest " << std::setw(3) << ++test_counter << ":  LLT factorization on " << symMatrix::name() << std::endl;

			auto A = symMatrix::get_random_matrix(sz, sz, 0.0, 1.0);

			// the rhs term of the linear system
			auto B = Matrix::get_random_matrix(sz, nrhs, 1.0, 0.5);

			for (index_t i = 1; i <= A.nrows(); ++i)
				A(i, i) += 100.0+i;

			Matrix X;
			symMatrix Acp(A);

			// executes the LLT factorization of Acp 
			timer.start();
			LLTSolver<symMatrix>::type  LLT(Acp);
			double t0 = timer.now();
			if (!LLT)
				std::cout << "Errore:" << LLT.state() << std::endl;
			else
			{
				// solves the linear system A X=B
				LLT.solve(X = B);
				double t1 = timer.now();
				if ( !LLT)
				{
					std::cout << "Errore\n";
					test_failure++;
				}
				else
				{
					// verify the correctness of the solution 
					Matrix tmp;
					MML3::Algorithm::product(A, X, tmp);
					tmp -= B;
					double err = MML3::Algorithm::norm2(tmp) / MML3::Algorithm::norm2(B);
					if (err > 1e-14)
					{
						std::cout << "Error\n";
						test_failure++;
					}
					else
						std::cout << "OK"
						<< "\n factorization time : " << t0
						<< "\n solution time      : " << (t1 - t0) << std::endl;
				
				}
			}

		}


		{

			typedef Mat<double>::General::Rectangular::ColMajor::type Matrix;
			typedef Mat<double>::Symmetric::Rectangular::ColMajor::type symMatrix;
			std::cout << "\nTest " << std::setw(3) << ++test_counter << ":  LLT factorization on " << symMatrix::name() << std::endl;

			auto A = symMatrix::get_random_matrix(sz, sz, 0.0, 1.0);

			// the rhs term of the linear system
			auto B = Matrix::get_random_matrix(sz, nrhs, 1.0, 0.5);

			for (index_t i = 1; i <= A.nrows(); ++i)
				A(i, i) += 100.0 + i;

			Matrix X;
			symMatrix Acp(A);

			// executes the LLT factorization of Acp 
			timer.start();
			LLTSolver<symMatrix>::type  LLT(Acp);
			double t0 = timer.now();
			if (!LLT)
				std::cout << "Errore:" << LLT.state() << std::endl;
			else
			{
				// solves the linear system A X=B
				LLT.solve(X = B);
				double t1 = timer.now();
				if (!LLT)
				{
					std::cout << "Errore\n";
					test_failure++;
				}
				else
				{
					// verify the correctness of the solution 
					Matrix tmp;
					MML3::Algorithm::product(A, X, tmp);
					tmp -= B;
					double err = MML3::Algorithm::norm2(tmp) / MML3::Algorithm::norm2(B);
					if (err > 1e-14)
					{
						std::cout << "Error\n";
						test_failure++;
					}
					else
						std::cout << "OK"
						<< "\n factorization time : " << t0
						<< "\n solution time      : " << (t1 - t0) << std::endl;

				}
			}

		}

		




		{

			typedef Mat<double>::Symmetric::LowerTria::ColMajor::type symMatrix;
			typedef Mat<double>::General::Rectangular::ColMajor::type Matrix;

			std::cout << "\nTest " << std::setw(3) << ++test_counter << ":  LLT factorization on " << symMatrix::name() << std::endl;


		auto A = symMatrix::get_random_matrix(sz, sz, 0.0, 0.0);

		for (index_t i = 1; i <= A.nrows(); ++i)
			A(i, i) = 100.0+i;


		//std::cout << "A packed:" << A << std::endl;

		// a copy of A
		auto Ac = A;
		// rhs terms
		auto B = Matrix::get_random_matrix(sz, nrhs, 0.0, 0.5);
		// solutions
		Matrix X(B), Bb(B);

		double t0,t1;
		//  solver for SYM matrice
		timer.start();
		LLTSolver<symMatrix>::type LLT(A);
		t0 = timer.now();
		
		if (LLT)
		{
			LLT.solve(X);
			t1 = timer.now();
			// hand made prod A*X
			for (index_t i = 1; i <= A.nrows(); ++i)
			{
				for (index_t k = 1; k <= X.ncols(); ++k)
				{

					double tmp = 0;
					for (int j = 1; j <=A.ncols(); ++j)
					{
						if ( j<=i)
							tmp += Ac(i, j)*X(j, k);
						else
							tmp += Ac(j, i)*X(j, k);
					}
					Bb(i, k) = tmp;
				}
			}


			Bb -= B;
			double err_sym = MML3::Algorithm::norm2(Bb) / MML3::Algorithm::norm2(B);
			
			if (err_sym > 1e-14)
			{
				std::cout << "Error\n";
				test_failure++;
			}
			else
			{
				std::cout << "OK"
					<< "\nfactor time : " << t0
					<< "\nsolver time : " << t1 - t0 << std::endl;
			}
		}
		else
		{
			std::cout << "Error\n";
			test_failure++;
		}

	}


		{

			
			typedef Mat<double>::Symmetric::UpperTria::ColMajor::type  symMatrix;
			typedef Mat<double>::General::Rectangular::ColMajor::type  Matrix;
			std::cout << "\nTest " << std::setw(3) << ++test_counter << ":  LLT factorization on " << symMatrix::name() << std::endl;
			auto A = symMatrix::get_random_matrix(sz, sz, 0.0, 5.0);

			for (index_t i = 1; i <= A.nrows(); ++i)
				A(i, i) = 100.0 + i;


			//std::cout << "A packed:" << A << std::endl;

			// a copy of A
			auto Ac = A;
			// rhs terms
			auto B = Matrix::get_random_matrix(sz, nrhs, 0.0, 0.5);
			// solutions
			Matrix X(B), Bb(B);

			double t0, t1, t2;
			//  solver for SYM matrice
			t0 = timer.start();
			LLTSolver<symMatrix>::type LLT(A);
			t1 = timer.now();

			if (LLT)
			{
				LLT.solve(X);
				t2 = timer.now();
				// hand made prod A*X
				for (int i = 1; i <= A.nrows(); ++i)
				{
					for (index_t k = 1; k <= X.ncols(); ++k)
					{

						double tmp = 0;
						for (int j = 1; j <= A.ncols(); ++j)
						{
							if (j >= i)
								tmp += Ac(i, j)*X(j, k);
							else
								tmp += Ac(j, i)*X(j, k);
						}
						Bb(i, k) = tmp;
					}
				}


				Bb -= B;
				double err_sym = MML3::Algorithm::norm2(Bb);
				if (err_sym > 1e-12)
					std::cout << "error on sym:" << err_sym << std::endl;
				if (err_sym > 1e-13)
				{
					std::cout << "Error\n";
					test_failure++;
				}
				else
				{
					std::cout << "OK"
						<< "\nfactor time : " << t1
						<< "\nsolver time : " << t2 - t1 << std::endl;

				}
			}
			else
			{
				std::cout << "Error\n";
				test_failure++;
			}



		}


	
		{

		

		typedef Mat<double>::Symmetric::LowerTria::RowMajor::type symMatrix;
		typedef Mat<double>::General::Rectangular::RowMajor::type  Matrix;
		std::cout << "\nTest " << std::setw(3) << ++test_counter << ":  LLT factorization on " << symMatrix::name() << std::endl;
		auto A = symMatrix::get_random_matrix(sz, sz, 0.0, 0.0);

		for (index_t i = 1; i <= A.nrows(); ++i)
			A(i, i) = 100.0 + i;


		//std::cout << "A packed:" << A << std::endl;

		// a copy of A
		auto Ac = A;
		// rhs terms
		auto B = Matrix::get_random_matrix(sz, nrhs, 0.0, 0.5);
		// solutions
		Matrix X(B), Bb(B);

		double t0, t1, t2;
		//  solver for SYM matrice
		t0 = timer.start();
		LLTSolver<symMatrix>::type LLT(A);
		t1 = timer.now();

		if (LLT)
		{
			LLT.solve(X);
			t2 = timer.now();
			// hand made prod A*X
			for (int i = 1; i <= A.nrows(); ++i)
			{
				for (index_t k = 1; k <= X.ncols(); ++k)
				{

					double tmp = 0;
					for (int j = 1; j <= A.ncols(); ++j)
					{
						if (j <= i)
							tmp += Ac(i, j)*X(j, k);
						else
							tmp += Ac(j, i)*X(j, k);
					}
					Bb(i, k) = tmp;
				}
			}


			Bb -= B;
			double err_sym = MML3::Algorithm::norm2(Bb);
			if (err_sym > 1e-12)
				std::cout << "error on sym:" << err_sym << std::endl;
			if (err_sym > 1e-13)
			{
				std::cout << "Error\n";
				test_failure++;
			}
			else
			{
				std::cout << "OK"
					<< "\nfactor time : " << t1
					<< "\nsolver time : " << t2 - t1 << std::endl;

			}
		}
		else
		{
			std::cout << "Error\n";
			test_failure++;
		}



	}


		{



			//typedef Mat<double>::Symmetric::UT symMatrix;
			typedef Mat<double>::Symmetric::UpperTria::type symMatrix;
			typedef Mat<double>::type Matrix;
			std::cout << "\nTest " << std::setw(3) << ++test_counter << ":  LLT factorization on " << symMatrix::name() << std::endl;
			auto A = symMatrix::get_random_matrix(sz, sz, 0.0, 0.0);

			for (index_t i = 1; i <= A.nrows(); ++i)
				A(i, i) = 100.0 + i;


			//std::cout << "A packed:" << A << std::endl;

			// a copy of A
			auto Ac = A;
			// rhs terms
			auto B = Matrix::get_random_matrix(sz, nrhs, 0.0, 0.5);
			// solutions
			Matrix X(B), Bb(B);

			double t0, t1, t2;
			//  solver for SYM matrice
			t0 = timer.start();
			LLTSolver<symMatrix>::type LLT(A);
			t1 = timer.now();

			if (LLT)
			{
				LLT.solve(X);
				t2 = timer.now();
				// hand made prod A*X
				for (int i = 1; i <= A.nrows(); ++i)
				{
					for (index_t k = 1; k <= X.ncols(); ++k)
					{

						double tmp = 0;
						for (int j = 1; j <= A.ncols(); ++j)
						{
							if (j <= i)
								tmp += Ac(i, j)*X(j, k);
							else
								tmp += Ac(j, i)*X(j, k);
						}
						Bb(i, k) = tmp;
					}
				}


				Bb -= B;
				double err_sym = MML3::Algorithm::norm2(Bb);
				if (err_sym > 1e-12)
					std::cout << "error on sym:" << err_sym << std::endl;
				if (err_sym > 1e-13)
				{
					std::cout << "Error\n";
					test_failure++;
				}
				else
				{
					std::cout << "OK"
						<< "\nfactor time : " << t1
						<< "\nsolver time : " << t2 - t1 << std::endl;

				}
			}
			else
			{
				std::cout << "Error\n";
				test_failure++;
			}



		}

	


		std::cout << "\n\nTest performed: " << test_counter << std::endl << "Test failed   : " << test_failure << std::endl;


		std::cout << "Elapsed time (sec): " << timer.seconds() << std::endl
			<< "Accuracy     (sec): " << std::setprecision(4) << timer.accuracy() << std::endl;
	}
	catch (std::exception& e)
	{

		std::cout << "Exception: " << e.what() << std::endl;
		return -1;
	}



return 0;
}

