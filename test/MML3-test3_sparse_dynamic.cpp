


#include"MML3-Timer.h"
#include"MML3-Matrix.h"
#include"MML3-dynamic_sparse_matrix.h"
#include"MML3-CSR3_Matrix.h"
#include"MML3-RangeT.h"
#include<iostream>

using namespace MML3;



	int main(){
		//typedef dynamic_sparse_row_Matrix<double, int, M_PROP::SYM, ordered_list_slist<int, double>> dsMatrix;
		typedef sparse_sym_matrix<double, int> dsMatrix;
		typedef CSR3_sym_matrix<double, int>		   ssMatrix;
		
		
		try
		{
			int sz = 4000;
			
			dsMatrix A = dsMatrix::get_random_matrix(sz, sz, 1. / sz, 10.0, 50.);

			ssMatrix sA(A.nrows(), A.ncols(), A.nonzeros());

			A.copy2CSR3(sA.nrows() + 1, sA.nonzeros(), sA.row_pos(), sA.column_index(), sA.column_value());
	//		std::cout << "sA: " << sA << std::endl;

			A.sum(1, 1, 1.0);
			A.sum(123, 50, 123);
			A.sum(50, 128, 128.0);
			A.put2CSR3(sA.nrows() + 1, sA.nonzeros(), sA.row_pos(), sA.column_index(), sA.column_value());
	//		std::cout << "sA: " << sA << std::endl;


			/*
			MML3::Timer timer;
			timer.start();
			// makes B a sz*sz sparse matrix with a mean of 40 elements per row
			auto B = dsMatrix::get_random_matrix(sz, sz, 40. / sz, 10.0, 50.);
			double t0 = timer.now();

			size_t nnz = B.nonzeros();
			std::cout	<< "size    =" << sz 
						<< "\tnnz=" << nnz << std::endl
						<< "sparsity=" << nnz / (sz*double(sz)) << std::endl
						<< "time    =" << t0 << std::endl
						<< "time/nnz=" << t0 / nnz << std::endl;
			
			int rv=B.fwrite("pippo.smat");
			B.destroy();
			rv = B.fread("pippo.smat");
			B.resize(1000, 1000);
			
			auto  k =  dsMatrix::inserter_matrix_t::get_random_matrix(12, 12, 10., 15.);
			auto  ks = dsMatrix::inserter_sym_matrix_t::get_random_matrix(12, 12, 10., 15.);
			dsMatrix::iset_t	ir = { 1, 3, 6, 8, 9,23,25,789,790,791,998,1000 };
			B.sum(ir, ir, k);
			std::cout << "B" << B << std::endl;
			B.sum(ir, ks);
			std::cout << "B" << B << std::endl;

			// iterating trough the matrix elements
			for (int r = 1; r <= B.nrows(); ++r)
			{
				// auto is a dsMatrix::iterator
				for (auto rit=B.row_begin(r); rit != B.row_end(r); ++rit)
				{
					//rit->idx   is the component column
					//rit->value is the component value
					std::cout << "[" << r << "," << rit->first <<  "]=" << rit->second << std::endl;
				}
			}
			
			*/
			



		}
		catch (std::exception& exc)
		{

			std::cout << "EXCEPTION : " << exc.what() << std::endl;
		}
		
		
		return 0;

		
	}

