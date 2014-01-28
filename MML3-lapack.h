#ifndef __MML3_CPP_LAPACK_H__
#define __MML3_CPP_LAPACK_H__
// MODIFY HERE TO INCLUDE YOUR OWN COMPLEX TYPE DEFINITIONS
#include<complex>
#include"MML3-config.h"


#if defined(MML3_LAPACK)
	// intel versione
	#if MML3_LAPACK==1
		#define MKL_Complex8		MML3::complex<float>
		#define MKL_Complex16		MML3::complex<double>
		#include "mkl_lapacke.h"
		static_assert(std::is_same<MML3::lapack_int_t, MKL_INT>::value, " MML3 Lapack int type differs from MKL int type");
	#elif MML3_LAPACK==2 // netlib version
		#define DLAPACK_COMPLEX_CPP
		#define lapack_complex_float  MML3::complex<float>
		#define lapack_complex_double MML3::complex<double>
		//lapack_complex_float lapack_make_complex_float(float re, float im);
		//lapack_complex_double lapack_make_complex_double(double re, double im);
		#include "lapacke.h"
		static_assert(std::is_same <MML3::lapack_int_t, lapack_int>::value, " MML3 Lapack int type differs from netlib lapack int type");
	#else
		#error   CANNOT COMPILE WITH THE SPECIFIED VALUE FOR  MML3_LAPACK= in  file MML3-config
	#endif

	
#endif



// ===========================================================================
// Class  BLAS1, BLAS2 , BLAS3  FUNCTOR CLASSES
// C++ wrapper to the cblas functions.
// ===========================================================================
namespace MML3{


		template<typename T>
		struct LAPACK
		{
			typedef	T					value_t;
			typedef lapack_int_t		int_t;
			
			// FACTORIZATION ROUTINES
			//  LU for dense matrices
			static int  getrf(int_t matrix_order, int_t m, int_t n, value_t* a, int_t lda, int_t* ipiv);
			//  LLT for dense symmetric matrices
			static int  potrf(int_t matrix_order, char uplo, int_t n, value_t* a, int_t lda);
			//  LLT for packed  symmetric matrices
			static int  pptrf(int_t matrix_order, char uplo, int_t n, value_t *ap);

			//  SOLUTION ROUTINES 
			//  LU  for dense matrices
			static int  getrs(int_t matrix_order, char trans, int_t n, int_t nrhs, const value_t* a, int_t lda, int_t* ipiv, value_t* b, int_t ldb);
			//  LLT for dense symmetric matrices
			static int  potrs(int_t matrix_order, char uplo, int_t n, int_t nrhs, const value_t* a, int_t lda, value_t* b, int_t ldb);
			//  LLT for packed  symmetric matrices
			static int  pptrs(int_t matrix_order, char uplo, int_t n, int_t nrhs, const value_t* ap, value_t* b, int_t ldb);
				   
			
		};


		///////////////////
		/// GETRF
		///////////////////
		template<>
		int LAPACK<float>::getrf(int matrix_order, int m, int n, value_t* a, int lda, int_t* ipiv){ return LAPACKE_sgetrf(matrix_order, m, n, a, lda, ipiv); }
		template<>
		int LAPACK<double>::getrf(int matrix_order, int m, int n, value_t* a, int lda, int_t* ipiv){ return LAPACKE_dgetrf(matrix_order, m, n, a, lda, ipiv); }
		template<>
		int LAPACK<complex<float>>::getrf(int matrix_order, int m, int n, value_t* a, int lda, int* ipiv){ return LAPACKE_cgetrf(matrix_order, m, n, a, lda, ipiv); }
		template<>
		int LAPACK<complex<double>>::getrf(int matrix_order, int m, int n, value_t* a, int lda, int* ipiv){ return LAPACKE_zgetrf(matrix_order, m, n, a, lda, ipiv); }

		///////////////////
		/// POTRF
		///////////////////
		template<>
		int LAPACK<float>::potrf(int matrix_order, char uplo, int n, value_t* a, int lda){ return LAPACKE_spotrf(matrix_order, uplo, n, a, lda); }
		template<>
		int LAPACK<double>::potrf(int matrix_order, char uplo, int n, value_t* a, int lda){ return LAPACKE_dpotrf(matrix_order, uplo, n, a, lda); }
		template<>
		int LAPACK<complex<float>>::potrf(int matrix_order, char uplo, int n, value_t* a, int lda){ return LAPACKE_cpotrf(matrix_order, uplo, n, a, lda); }
		template<>
		int LAPACK<complex<double>>::potrf(int matrix_order, char uplo, int n, value_t* a, int lda){ return LAPACKE_zpotrf(matrix_order, uplo, n, a, lda); }

		///////////////////
		/// GETRS
		///////////////////
		template<>
		int LAPACK<float>::getrs(int matrix_order, char trans, int n, int nrhs, const value_t* a, int lda, int* ipiv, value_t* b, int ldb){
			return LAPACKE_sgetrs(matrix_order, trans, n, nrhs, a, lda, ipiv, b, ldb);
		}
		template<>
		int LAPACK<double>::getrs(int matrix_order, char trans, int n, int nrhs, const value_t* a, int lda, int* ipiv, value_t* b, int ldb){
			return LAPACKE_dgetrs(matrix_order, trans, n, nrhs, a, lda, ipiv, b, ldb);
		}

		template<>
		int LAPACK<complex<float>>::getrs(int matrix_order, char trans, int n, int nrhs, const value_t* a, int lda, int* ipiv, value_t* b, int ldb){
			return LAPACKE_cgetrs(matrix_order, trans, n, nrhs, a, lda, ipiv, b, ldb);
		}
		template<>
		int LAPACK<complex<double>>::getrs(int matrix_order, char trans, int n, int nrhs, const value_t* a, int lda, int* ipiv, value_t* b, int ldb){
			return LAPACKE_zgetrs(matrix_order, trans, n, nrhs, a, lda, ipiv, b, ldb);
		}



		///////////////////
		/// POTRS
		///////////////////
		template<>
		int LAPACK<float>::potrs(int matrix_order, char uplo, int n, int nrhs, const value_t* a, int lda, float* b, int ldb){
			return LAPACKE_spotrs(matrix_order, uplo, n, nrhs, a, lda, b, ldb);
		}
		template<>
		int LAPACK<double> ::potrs(int matrix_order, char uplo, int n, int nrhs, const value_t* a, int lda, double* b, int ldb){
			return LAPACKE_dpotrs(matrix_order, uplo, n, nrhs, a, lda, b, ldb);
		}

		template<>
		int LAPACK<complex<float>> ::potrs(int matrix_order, char uplo, int n, int nrhs, const value_t* a, int lda, complex<float>* b, int ldb){
			return LAPACKE_cpotrs(matrix_order, uplo, n, nrhs, a, lda, b, ldb);
		}
		template<>
		int LAPACK<complex<double>>::potrs(int matrix_order, char uplo, int n, int nrhs, const value_t* a, int lda, complex<double>* b, int ldb){
			return LAPACKE_zpotrs(matrix_order, uplo, n, nrhs, a, lda, b, ldb);
		}

		///////////////////
		/// PPTRF
		///////////////////
		template<>
		int LAPACK<float>::pptrf(int matrix_order, char uplo, int n, float* ap){
			return LAPACKE_spptrf(matrix_order, uplo, n, ap);
		}
		template<>
		int LAPACK<double>::pptrf(int matrix_order, char uplo, int n, double* ap){
			return LAPACKE_dpptrf(matrix_order, uplo, n, ap);
		}
		template<>
		int LAPACK<complex<float>>::pptrf(int matrix_order, char uplo, int n, complex<float>* ap){
			return LAPACKE_cpptrf(matrix_order, uplo, n, ap);
		}
		template<>
		int LAPACK<complex<double>>::pptrf(int matrix_order, char uplo, int n, complex<double>* ap){
			return LAPACKE_zpptrf(matrix_order, uplo, n, ap);
		}

		///////////////////
		/// PPTRS
		///////////////////
		template<>
		int LAPACK<float>::pptrs(int matrix_order, char uplo, int n,	int nrhs, const float* ap, float* b,int ldb) {
			return LAPACKE_spptrs(matrix_order, uplo, n,nrhs, ap, b,ldb);}

		template<>
		int LAPACK<double>::pptrs(int matrix_order, char uplo, int n, int nrhs, const double* ap, double* b, int ldb) {
			return LAPACKE_dpptrs(matrix_order, uplo, n, nrhs, ap, b, ldb);
		}

		template<>
		int LAPACK< complex<float>>::pptrs(int matrix_order, char uplo, int n, int nrhs, const complex<float>* ap, complex<float>* b, int ldb) {
			return LAPACKE_cpptrs(matrix_order, uplo, n, nrhs, ap, b, ldb);
		}
		template<>
		int LAPACK< complex<double>>::pptrs(int matrix_order, char uplo, int n, int nrhs, const complex<double>* ap, complex<double>* b, int ldb) {
			return LAPACKE_zpptrs(matrix_order, uplo, n, nrhs, ap, b, ldb);
		}

		

}
#endif
