

#include"MML3-Timer.h"
#include"MML3-Matrix.h"
#include"MML3-MatrixAlgorithms.h"
#include"MML3-dynamic_CSR_Matrix.h"
#include"MML3-static_CSR_Matrix.h"
#include"MML3-RangeT.h"
#include<iostream>

using namespace MML3;


	int main(){

		try
		{
			
			int sz = 20;
			double sparsity = 2. / sz;
			typedef dynamic_sparse_sym_matrix<double>  dynaMatrix;
			dynaMatrix    A=dynaMatrix::get_random_matrix(sz, sz, sparsity, 10.0, 50.);

			int rv = A.fwrite("pippo.smat");
			if (rv!=0)
				throw std::runtime_error("Errore scrittura dyn");

			typedef static_sparse_sym_matrix<double>  statSymMatrix;
			
			statSymMatrix B;
			rv=B.fread("pippo.smat");
			if (rv!=0)
				throw std::runtime_error("Errore lettura stat");

			// confronto le due ,matrici usando gli iteratori sugli elementi delle righe
			double tmp = 0.;
			for (auto r : Range<>(1, sz))
				for (auto it = A.row_begin(r), end = A.row_end(r); it != end; ++it)
					tmp += fabs(*B.get_p(r, it.idx()) - it.value());
			
			std::cout << "diff=" << tmp << std::endl;
			if (tmp > 1E-15)
				throw std::runtime_error("Errore differenza");

			
			// anche le statiche hanno gli iteratori sugli elementi di una riga
			for (auto r : Range<>(1, sz))
				for (auto it = B.row_begin(r), end = B.row_end(r); it != end; ++it)
					std::cout << r << "\t" << it.idx() << "\t" << it.value() << std::endl;
			{

					typedef dynamic_sparse_matrix<double>		dyGeMatrix;
					typedef dynamic_sparse_sym_matrix<double>	dySymMatrix;

					dyGeMatrix A(6, 6);
					A.put(1, 1, 11.0);
					A.put(1, 6, 16.0); A.put(6, 1, 16.0);
					A.put(2, 3, 23.0); A.put(3, 2, 23.0);
					A.put(2, 5, 25.0); A.put(5, 2, 25.0);
					A.put(3, 3, 33.0);
					A.put(4, 5, 45.0); A.put( 5,4, 45.0);
					A.put(5, 6, 56.0); A.put(6, 5, 56.0);
					std::cout << A << std::endl;

					dySymMatrix As(6, 6);
					As.put(1, 1, 11.0);
					As.put(1, 6, 16.0); As.put(6, 1, 16.0);
					As.put(2, 3, 23.0); As.put(3, 2, 23.0);
					As.put(2, 5, 25.0); As.put(5, 2, 25.0);
					As.put(3, 3, 33.0);
					As.put(4, 5, 45.0); As.put(5, 4, 45.0);
					As.put(5, 6, 56.0); As.put(6, 5, 56.0);
					std::cout << As << std::endl;

					int rv=A.fwrite("A.smat");
					rv = As.fwrite("As.smat");


					typedef static_sparse_matrix<double>		stGeMatrix;
					typedef static_sparse_sym_matrix<double>	stSymMatrix;
					stGeMatrix B;
					stSymMatrix Bs;
					rv=B.fread("A.smat");
					rv = Bs.fread("As.smat");
					std::cout << B << std::endl;
					std::cout << Bs << std::endl;

					Mat<double>::type   b = Mat<double>::type::get_random_matrix(6, 1, 1.0, 0.5), c(6), c1(6);
					b.fill(1.0);

					B.product(b.begin(), b.size(), c.begin(), c.size());
					Bs.product(b.begin(), b.size(), c1.begin(), c1.size());
					double diff = MML3::Algorithm::norm2(c - c1);
					if (diff > 1e-15)
					{
						std::cerr << "errore in product" << std::endl;
					}
			}

			{

				int sz = 6000;
				double sparsity = 50. / sz;
				typedef dynamic_sparse_sym_matrix<double>  dynaMatrix;
				std::cout	<< "building sparse symmetrix matrix of size " << sz << "\t sparsity: " << sparsity 
							<< "\t nnz : " << (sparsity*sz*sz)/1E6 << " M" << std::endl;
				dynaMatrix    A = dynaMatrix::get_random_matrix(sz, sz, sparsity, 10.0, 50.);
				std::cout << "writing  on file " << std::endl;
				int rv = A.fwrite("pippo.smat");
				if (rv != 0)
					throw std::runtime_error("Errore scrittura dyn");

				typedef static_sparse_sym_matrix<double>  statSymMatrix;
				statSymMatrix B;
				std::cout << "reading  from file " << std::endl;
				rv = B.fread("pippo.smat");
				if (rv != 0)
					throw std::runtime_error("Errore lettura stat");

				Mat<double>::type  b(sz), c(sz);
				b.fill(3.14);

				std::cout << "matrix vector product " << std::endl;
				MML3::Timer timer;
				timer.start();
				B.product(b.begin(), b.size(), c.begin(), c.size());
				double t0 = timer.now();
				std::cout << "size: " << sz << "\tsparsity: " << double(B.nonzeros()) / (double(sz)*sz)
					<< "\ttime: " << t0 << std::endl;
			}

		}
		catch (std::exception& exc)
		{
			std::cout << "EXCEPTION : " << exc.what() << std::endl;
		}
		
		
		return 0;
	}

