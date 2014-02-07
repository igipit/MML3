
#include"MML3-Vector.h"
#include"MML3-Matrix.h"
#include"mml3-sparseeigensolver.h"
#include"mkl.h"
#include"mkl_dss.h"
#include"mkl_pardiso.h"
#include"mkl_spblas.h"
#include<map>
#include<array>


namespace MML3
{

	// solutore per la determinazione dei primi modi propri generalizzati (l_i, Z_i) per matrici CSR3  simmetriche 
	// A Z_i = li B Z_i
	int SparseEigenSolver(	static_sparse_CSR_Matrix<double, std::int32_t, M_PROP::SYM>& A,
							static_sparse_CSR_Matrix<double, std::int32_t, M_PROP::SYM>& B,
							Matrix<double, M_PROP::GE, M_SHAPE::RE, M_ORD::COL>& eigen_vector,
							Vector<double>& eigen_value,
							Vector<double>& eigen_error,
							double eigenv_max,
							std::array<int,32>& param,
							std::ostream& log)
{

		typedef MML3_INT_TYPE  int_t;
		int len=100;
		char buf[100];
		//MKLGetVersionString(buf, len);
		 mkl_get_version_string(buf, len);

		log	<< "Solutore modale per matrici sparse basato su:\n"
			<< "\tMKL-FEAST: " << std::string(buf,len);
		

		log << "\nCaratteristiche del sistema:"
			<< "\n\t    N      =" << A.nrows()
			<< "\n\t    NNZ(A) =" << A.nonzeros()
			<< "\n\t    NNZ(B)="  << B.nonzeros() << std::endl;

		if (A.nrows() != B.nrows())
		{
			log << "LinearEigenSolver: dimensioni di A e B incompatibili \n";
			return -1;
		}


		Vector<int_t> fpm(128);
		feastinit(fpm.begin());

		fpm(1) = 1; // manda in output messAGGI a runtime 
		fpm(2) = 8; // numero di punti di quadratura (default=)8
		fpm(3) = 12; // tolleranza 10^-fpm(3)
		fpm(64) = 0; // use PARDISO

		int_t	n				= A.nrows();
		double	epsout			= 0.;
		int_t	loop			= 0;
		double	emin			= 0.;
		int_t	subspace_dim	= 3; // m0 in feast routine
		
		int_t					founded_eig = 0; // feast m
		int_t					info = 0;
		const char uplo = 'U';
				
		dfeast_scsrgv(&uplo, &n, 
			A.column_value(), A.row_pos(), A.column_index(), 
			B.column_value(), B.row_pos(), B.column_index(), 
			fpm.begin(), &epsout, &loop, &emin, &eigenv_max,
			&subspace_dim, eigen_value.begin(), eigen_vector.begin(), &founded_eig, eigen_error.begin(), 
			&info);

		// successo
		if (info == 0)
		{


		}

		

	return 0;
		
}


} // end namespace XSA