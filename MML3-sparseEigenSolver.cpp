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
	
	int SparseEigenSolver(	
		CSR3_sym_matrix<double, std::int32_t>& A,
		CSR3_sym_matrix<double, std::int32_t>& B,
		Matrix<double, M_PROP::GE, M_SHAPE::RE, M_ORD::COL>& Z,
		Vector<double>& L,
		Vector<double>& eigen_error,
		double eigenv_max,
		int    toll_exp,
		int    max_iter,
		std::ostream& log)
{

		typedef MML3_INT_TYPE  int_t;
		int len=100;
		char buf[100];
		//MKLGetVersionString(buf, len);
		 mkl_get_version_string(buf, len);

		log	<< "\tSolutore modale per matrici sparse basato su:\n"
			<< "\tMKL-FEAST: " << std::string(buf,len);
		

		log << "\n\tCaratteristiche del sistema:"
			<< "\n\t\tN      : " << A.nrows()
			<< "\n\t\tNNZ(A) : " << A.nonzeros()
			<< "\n\t\tNNZ(B) : " << B.nonzeros() << std::endl;

		if (A.nrows() != B.nrows())
		{
			log << "LinearEigenSolver: dimensioni di A e B incompatibili \n";
			return -1;
		}



		


		Vector<int_t> fpm(128);
		feastinit(fpm.data());

		fpm(1) = 1; // manda in output messAGGI a runtime 
		fpm(2) = 8; // numero di punti di quadratura (default=)8
		if (toll_exp)
			fpm(3) =toll_exp; // tolleranza 10^-fpm(3)
		if (max_iter)
			fpm(4) = max_iter;
		fpm(64) = 0; // use PARDISO

		int_t	n				= A.nrows();
		double	epsout			= 0.;
		int_t	loop			= 0;
		double	emin			= 0.;
		int_t	subspace_dim	= std::min(Z.ncols(),n); // m0 in feast routine
		
		int_t					founded_eig = 0; // feast m
		int_t					info = 0;
		const char uplo = 'U';


		log << "\n\tParmetri algoritmici:"
			<< "\n\t\tintervallo autovalori     : [" << emin << "," << eigenv_max << "]"
			<< "\n\t\tdimensione sottospazio    : "		<< subspace_dim 
			<< "\n\t\ttolleranza                : 10E-" << fpm(3) 
			<< "\n\t\titerazioni massime        : "		<< fpm(4)
			<< std::endl;
		log << "\n ................ Inizio messaggi solutore .................\n\n";
		dfeast_scsrgv(&uplo, &n, 
			A.column_value(), A.row_pos(), A.column_index(), 
			B.column_value(), B.row_pos(), B.column_index(), 
			fpm.data(), &epsout, &loop, &emin, &eigenv_max,
			&subspace_dim, L.data(), Z.data(), &founded_eig, eigen_error.data(), 
			&info);
		log << "\n ................ Fine messaggi solutore ...................\n\n";


		if (info != 0)
			return info;

		Vector<double> tmp_L(L.data(), founded_eig);
		L.swap(tmp_L);

	return info;
		
}


} // end namespace XSA