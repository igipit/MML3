
#include"MML3-Vector.h"
#include"MML3-Matrix.h"
#include"MML3-sparseLinearSolver.h"
#include"mkl.h"
#include"mkl_dss.h"
#include"mkl_pardiso.h"
#include"mkl_spblas.h"


namespace MML3
{

	int LinearSparseSolver( static_sparse_CSR_Matrix<double, std::int32_t, M_PROP::SYM>& A,
							Matrix<double, M_PROP::GE, M_SHAPE::RE, M_ORD::COL>& B,
							Matrix<double, M_PROP::GE, M_SHAPE::RE, M_ORD::COL>& X,
							std::ostream& log,
							int threads)
{

		typedef MML3_INT_TYPE  int_t;
		
		int len=100;
		char buf[100];
		//MKLGetVersionString(buf, len);
		 mkl_get_version_string(buf, len);

		log	<< "Solutore lineare per matrici sparse basato su:\n"
			<< "\tMKL-PARDISO: " << std::string(buf,len);
		
		

		log << "\nCaratteristiche del sistema:"
			<< "\n\t    N   =" << A.nrows()
			<< "\n\t    NNZ =" << A.nonzeros()
			<< "\n\t    NRHS=" << B.ncols() << std::endl;
			
	// pardiso parameters
	MML3::Vector<int_t>	pardiso_pt(64);
	pardiso_pt=0;
	int_t				pardiso_maxfct = 1;	// una sola matrice da fattorizzare
	int_t				pardiso_mnum = 1;	// la matrice numero 1
	int_t				pardiso_mtype = 2;	// matrici simmetriche definite positive
	int_t				pardiso_phase = 0;
	int_t				pardiso_n = A.nrows();
	double*				pardiso_a=A.column_value();
	int_t*				pardiso_ia = A.row_pos();
	int_t*				pardiso_ja = A.column_index();
	MML3::Vector<int_t>	pardiso_perm(pardiso_n);
	int_t				pardiso_nrhs = 0;
	MML3::Vector<int_t>	pardiso_iparam(64); 
	pardiso_iparam = 0;
	int_t				pardiso_msglvl = 0;	// 1=messaggi a run time, 0=no messaggi
	double*				pardiso_b=0;
	double*				pardiso_x=0;
	int_t				pardiso_error=0;



	
	pardiso_iparam(1)  = 1; /* No solver default*/
	pardiso_iparam(2)  = 2; /* Fill-in reordering from METIS, 1 for AMD */
	pardiso_iparam(3)  = threads; /* Numbers of thread processors */
	pardiso_iparam(4)  = 0; /* No iterative-direct algorithm */
	pardiso_iparam(5)  = 0; /* No user fill-in reducing permutation */
	pardiso_iparam(6)  = 0; /* Write solution into x */
	pardiso_iparam(7)  = 16; /* Default logical fortran unit number for output */
	pardiso_iparam(8)  = 0; /* Max numbers of iterative refinement steps */
	pardiso_iparam(9)  = 0; /* Not in use*/
	pardiso_iparam(10) = 13; /* Perturb the pivot elements with 1E-13 */
	pardiso_iparam(11) = 0; /* Use nonsymmetric permutation and scaling MPS */
	pardiso_iparam(12) = 0; /* Not in use*/
	pardiso_iparam(13) = 0; /* Not in use*/
	pardiso_iparam(14) = 0; /* Output: Number of perturbed pivots */
	pardiso_iparam(15) = 0; /* Not in use*/
	pardiso_iparam(16) = 0; /* Not in use*/
	pardiso_iparam(17) = 0; /* Not in use*/
	pardiso_iparam(18) = -1; /* Output: Number of nonzeros in the factor LU */
	pardiso_iparam(19) = -1; /* Output: Mflops for LU factorization */
	pardiso_iparam(20) = 0;  /* Output: Numbers of CG Iterations */




	
	// --------------------------------------------------------------------
	// Reordering and Symbolic Factorization. This step also allocates 
	// all memory that is necessary for the factorization. 
	// --------------------------------------------------------------------
	pardiso_phase = 11;
	PARDISO(	pardiso_pt.begin(), &pardiso_maxfct, &pardiso_mnum, &pardiso_mtype, &pardiso_phase, 
				&pardiso_n, pardiso_a, pardiso_ia, pardiso_ja, pardiso_perm.begin(), &pardiso_nrhs, 
				pardiso_iparam.begin(), &pardiso_msglvl, 0, 0, &pardiso_error);

	if (pardiso_error != 0) 
	{
		log << "riordino e fattorizzazione simbolica falliti, codice: " << pardiso_error << "\n";
		return -1;
	}
	log <<	"\tRiordino e fattorizzazione simbolica completata" 
		<<	"\n\t\tFactor NNZ=" <<  pardiso_iparam(18) 
		<<	"\n\t\tMFLOPS    ="	<<  pardiso_iparam(19);


	// --------------------------------------------------------------------
	// Numerical factorization.
	// --------------------------------------------------------------------
	pardiso_phase = 22;
	PARDISO(	pardiso_pt.begin(), &pardiso_maxfct, &pardiso_mnum, &pardiso_mtype, &pardiso_phase, 
				&pardiso_n, pardiso_a, pardiso_ia, pardiso_ja, pardiso_perm.begin(), &pardiso_nrhs, 
				pardiso_iparam.begin(), &pardiso_msglvl, 0, 0, &pardiso_error);

	if (pardiso_error != 0) 
	{
		log << " fattorizzazione numerica fallita, codice pardiso=" << pardiso_error << "\n";
		return -1;
	}
		

	// --------------------------------------------------------------------*/
	//  Back substitution and iterative refinement. */
	// --------------------------------------------------------------------*/
	pardiso_phase = 33;
	pardiso_nrhs=B.nrows();
	pardiso_iparam(6)  = 0;  // Write solution into x. Attenzione, anche se lo pongo ad 1 in modo che la soluzione sovrascriva F, X viene comunque usato
	
	
	PARDISO(	pardiso_pt.begin(), &pardiso_maxfct, &pardiso_mnum, &pardiso_mtype, &pardiso_phase, 
				&pardiso_n, pardiso_a, pardiso_ia, pardiso_ja, pardiso_perm.begin(), &pardiso_nrhs, 
				pardiso_iparam.begin(), &pardiso_msglvl, B.begin(), X.begin(), &pardiso_error);

	if (pardiso_error != 0) 
	{
		log << "Errore nella fase  di sostituzione all'indietro e rifinitura iteritava, codice: " << pardiso_error << "\n";
		return -1;
	}
	

	// --------------------------------------------------------------------
	// Termination and release of memory. 
	// --------------------------------------------------------------------
	pardiso_phase = -1; /* Release internal memory. */
	PARDISO(	pardiso_pt.begin(), &pardiso_maxfct, &pardiso_mnum, &pardiso_mtype, &pardiso_phase, 
				&pardiso_n, pardiso_a, pardiso_ia, pardiso_ja, pardiso_perm.begin(), &pardiso_nrhs, 
				pardiso_iparam.begin(), &pardiso_msglvl, 0, 0, &pardiso_error);


	return 1;
		
}


} // end namespace XSA