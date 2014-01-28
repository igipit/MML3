

#include"MML3-Timer.h"
#include"MML3-Matrix.h"
#include"MML3-dynamic_CSR_Matrix.h"
#include"MML3-RangeT.h"
#include<iostream>

using namespace MML3;



	int main(){

		try
		{
			int sz = 1000;
					
			typedef dynamic_sparse_sym_matrix<double>  dsMatrix;
			dsMatrix A(sz,sz);

			// inserting some components
			A.put(1, 1, 1.0);
			A.put(50, 123, 123);
			A.put(50, 128, 128.0);

			std::cout << "A: " << A << std::endl;;

			MML3::Timer timer;
			timer.start();
			// makes B a sz*sz sparse matrix with a mena of 20 elements per row
			auto B = dsMatrix::get_random_matrix(sz, sz, 20. / sz, 10.0, 50.);
			double t0 = timer.now();

			size_t nnz = B.nonzeros();
			std::cout	<< "size    =" << sz 
						<< "\tnnz=" << nnz << std::endl
						<< "sparsity=" << nnz / (sz*double(sz)) << std::endl
						<< "time    =" << t0 
						<< "time/nnz=" << t0 / nnz << std::endl;

			// iterate trough columns of row 2
			dsMatrix::const_col_iterator it		= B.row_begin(2);
			dsMatrix::const_col_iterator end	= B.row_end(2);
			while (it != end)
				++it;
			int rv=B.fwrite("pippo.smat");
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
					std::cout << "[" << r << "," << rit->idx <<  "]=" << rit->value << std::endl;
				}
			}


			




















		}
		catch (std::exception& exc)
		{

			std::cout << "EXCEPTION : " << exc.what() << std::endl;
		}
		
		
		return 0;
	}

