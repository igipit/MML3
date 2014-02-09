
#include"MML3-config.h"
#if defined(MML3_LAPACK)
	#include"MML3-lapack.h"
#else
#error MML3::Solver relies upon Lapacke
#endif


#pragma once
#include"MML3-Matrix.h"
#include<valarray>
#include<memory>

namespace MML3{
	namespace Solver{

		// Gauss factorization solver
		template<typename T, typename MP, typename MS, typename MO>
		class LU;
		// Cholesky factorization solver
		template<typename T, typename MP, typename MS, typename MO>
		class LLT;

		template<typename MAT>
		struct LUSolver
		{
			typedef typename MAT::value_t    MT;
			typedef typename MAT::property_t MP;
			typedef typename MAT::shape_t	MS;
			typedef typename MAT::order_t    MO;
			typedef LU<MT, MP, MS, MO> type;
		};

		template<typename MAT>
		struct LLTSolver
		{
			typedef typename MAT::value_t    MT;
			typedef typename MAT::property_t MP;
			typedef typename MAT::shape_t	MS;
			typedef typename MAT::order_t    MO;
			typedef LLT<MT, MP, MS, MO> type;
		};






		// errors codes in the range (-9, -1) are reserved for lapack errors : -i means that the i'th parameter passed to lapack routine is wrong
		// error codes in the range (1,...)  in the factorization routine indicates that the i'th equation cannot be factorized
		// 0 return codes means that its all right
		enum State { OK = 0, NO_INIT = -100, BAD_SZ = -10, BAD_ON_SOLVE = -11, BAD_MTYPE };

		template<typename T, typename MP, typename MS, typename MO>
		class LU
		{
		public:
			typedef Matrix<T,MP,MS,MO>		matrix_t;
			typedef T 						value_t;
			typedef	lapack_int_t			int_t;
			static_assert(matrix_t::IS::GE && matrix_t::IS::RE, " only General Rectangular matrices are supported by LU solver");


		public:
			LU(matrix_t& A){ factorize(A); }
			// factorizes the matrix A and mantains internally a reference to A and a pivot vector
			~LU(){pA_ = nullptr;}


			bool	factorize(matrix_t& A);
			bool	solve(matrix_t& X, bool trA = false);
			int		state()const{return state_;}
			bool	good()const{ return state_ == State::OK; }

			operator bool()const{ return good();}
		private:
			int state_ = NO_INIT;
			std::valarray<int_t> pivot_;
			matrix_t*   pA_ = nullptr;
		};





		//  LLT solver for symmetric positive definite  matrices
		template<typename T, typename MP, typename MS, typename MO>
		class LLT
		{
			typedef Matrix<T, MP, MS, MO>	matrix_t;
			typedef T 						value_t;
			typedef	lapack_int_t			int_t;


			//static_assert(MAT::IS::RE || (MAT::IS::Packed  && MAT::IS::ColMajor),
			//	"Packed Symmetric Row major matrices are not fully supported by Lapacke ?pptrf\s function (version 11.1)");
			// but i provided a workaround

		public:
			LLT(matrix_t& A, bool L = true) { factorize(A, L); }
			~LLT(){ pA_ = nullptr; }

			int  state()const{ return state_; }
			bool good()const{ return state_ == 0;}
			operator bool()const{ return good(); }

			bool factorize(matrix_t& A, bool L = true);

			template<typename BMAT>
			bool solve(BMAT& X);
		private:
			int state_ = NO_INIT;
			matrix_t*  pA_=nullptr;
			char uplo_;

		};










		//////////////////////////////////////////////
		// IMPLEMENTATION


		template<typename T, typename MP, typename MS, typename MO>
		bool LU<T,MP,MS,MO>::factorize(matrix_t& A)
		{
			pA_ = &A;
			pivot_.resize(A.nrows());
			state_=LAPACK<value_t>::getrf(A.CblasOrder(), A.nrows(), A.ncols(), A.data(), A.leading_dim(), &pivot_[0]);
			return good();
		}




		template<typename T, typename MP, typename MS, typename MO>
		bool LU<T, MP, MS, MO>::solve(matrix_t& X, bool trA)
		{
			if (state_ != 0)
				state_ = State::BAD_ON_SOLVE;
			else if (X.nrows() != pA_->nrows())
				state_ = State::BAD_SZ;
			else
			{
				char tr('N');
				if (trA)
					tr = 'S';
				state_ = LAPACK<value_t>::getrs(pA_->CblasOrder(), tr, pA_->nrows(), X.ncols(), pA_->data(), pA_->leading_dim(), &pivot_[0], X.data(), X.leading_dim());
			}
			return good();
		}



		template<typename T, typename MP, typename MS, typename MO>
		bool LLT<T,MP,MS,MO>::factorize(matrix_t& A, bool L )
		{
			uplo_ = L ? 'L' : 'U';
			pA_ = &A;
			// for rectangular matrices there is no problem for lapacke potrf
			if (std::is_same<MS,M_SHAPE::RE>::value)
			{
				// if is a rectangular symmetric the input paramater L is ignored
				if (std::is_same<MP,M_PROP::SYM>::value)
					uplo_ = (A.CblasSymUpLo() == Option::UpLo::Lower) ? 'L' : 'U';
				state_ = LAPACK<value_t >::potrf(A.CblasOrder(), uplo_, A.nrows(), A.data(), A.leading_dim());
			}
			// while for packed matrices it seems that intel mkl pptrf does not handle properly row major formats, so ...
			else
			{
				// Col major formats are treated correctly
				if (std::is_same<MO,M_ORD::COL>::value)
				{
					uplo_ = std::is_same<MS,M_SHAPE::LT>::value?'L':'U';
					state_ = LAPACK<value_t >::pptrf(A.CblasOrder(), uplo_, A.nrows(), A.begin());
				}
				// Row major formats require some workaround
				else
				{
					// I treat Row major  packed LT mats as they where ColMajor packed UT because the underlying linear vector is the same
					if (std::is_same<MS,M_SHAPE::LT>::value)
					{
						uplo_ = 'L';
						state_ = LAPACK<value_t >::pptrf(Option::ColMajor, 'U', A.nrows(), A.begin());
					}
					// Row major packed UT have the same array structures of Col major LT
					else
					{
						uplo_ = 'U';
						state_ = LAPACK<value_t >::pptrf(Option::ColMajor, 'L', A.nrows(), A.begin());
					}
				}
			}
			return good();
		}

		template<typename T, typename MP, typename MS, typename MO>
		template<typename BMAT>
		bool LLT<T, MP, MS, MO>::solve(BMAT& X)
		{
			if (state_ != 0)
				state_= State::BAD_ON_SOLVE;
			else if (X.nrows() != pA_->nrows())
				state_=State::BAD_SZ;
			else
			{
				if (std::is_same<MS,M_SHAPE::RE>::value)
					state_ = LAPACK<value_t>::potrs(pA_->CblasOrder(), uplo_, pA_->nrows(), X.ncols(), pA_->data(), pA_->leading_dim(), X.begin(), X.leading_dim());
				else
				{
					if (std::is_same<MO,M_ORD::COL>::value)
						state_ = LAPACK<value_t>::pptrs(pA_->CblasOrder(), uplo_, pA_->nrows(), X.ncols(), pA_->data(), X.begin(), X.leading_dim());
					else
					{
						char tmp_uplo = std::is_same<M_SHAPE::LT,MS>::value?'U':'L';

						// i make a transposed copy of the rhs
						auto B = transpose(X);
						state_ = LAPACK<value_t>::pptrs(Option::ColMajor, tmp_uplo, pA_->nrows(), X.ncols(), pA_->data(), B.begin(), B.ncols());
						X = transpose(B);
					}
				}
			}
			return good();
		}


		/*

		//  LLT solver specialized for Packed symmetric  matrices
		template<typename T>
		class LLT<Matrix<T, MAT_PROP::SYM>>
		{
			typedef Matrix<T, MAT_PROP::SYM>			sym_matrix_t;
			typedef typename sym_matrix_t::property_t		property_t;
			typedef typename sym_matrix_t::value_t 			value_t;
			typedef typename LAPACK<value_t>::int_t		int_t;
			typedef Matrix<T, MAT_PROP::GE>				rhs_matrix_t;
		public:
			LLT() = default;
			int factorize(sym_matrix_t& A)
			{
				pA_ = &A;
				// !!! con ROW_MAJOR ed 'L' non funziona e pertanto, in maniera equivalente, chiamo la routine
				// con COL_MAJOR e 'U'
				int info = LAPACK<T>::pptrf(MML3_LAPACK_COL_MAJOR, 'U', A.nrows(), A.begin());
				return info;
			}
			int solve(rhs_matrix_t& X)
			{
				if (X.nrows() != pA_->nrows())
					return std::numeric_limits<int>::min();
				// pare che pptrs non voglia saperne di funzionare in modalità ROW_MAJOR, pertanto eseguo la trasposta del termine noto
				X.set_transpose();
				T* px = X.begin();
				int info = LAPACK<T>::pptrs(MML3_LAPACK_COL_MAJOR, 'U', pA_->nrows(), X.nrows(), pA_->data(), px, X.leading_dimension());
				X.set_transpose();
				return info;
			}
		private:
			sym_matrix_t*  pA_{};
		};

		// an LLT maker
		template< typename T, typename MP>
		auto  make_LLT(const Matrix<T, MP>&)->LLT< Matrix<T, MP>>
		{
			return LLT<Matrix<T, MP>>();
		};
		*/



	} // end namespace
} // end namespace
