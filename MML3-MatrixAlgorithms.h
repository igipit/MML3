#pragma once
#include"MML3-config.h"
#include"MML3-Matrix.h"
#include"MML3-Vector.h"
#include"MML3-blas.h"
#include<limits>
#include<initializer_list>



namespace MML3
{
	namespace Algorithm
	{


		/////////////////////////
		// ALGORITHMS HEADERS  //
		/////////////////////////

		// y=y+a*x  ( BLAS: ?AXPY)
		template<typename T, typename MP, typename MS, typename MO>
		Matrix<T, MP, MS, MO>&		axpy(T a, const Matrix<T, MP, MS, MO>& x, Matrix<T, MP, MS, MO>& y);

		template<typename T>
		Vector<T>&					axpy(T a, const Vector<T>& x, Vector<T>& y);


		// x=x*a     (BLAS: ?SCAL)
		template<typename T, typename MP, typename MS, typename MO>
		Matrix<T, MP, MS, MO>&		scal(Matrix<T, MP, MS, MO>& x, T a);
		template<typename T>
		Vector<T>&					scal(Vector<T>& x, T a);


		// scalar product  ps(A)= \sum_{i,j} A_ij* A_ij
		template<typename T, typename MP, typename MS, typename MO>
		T							dot(const Matrix<T, MP, MS, MO>& x, const Matrix<T, MP, MS, MO>& y);

		template<typename T>
		T							dot(const Vector<T>& x, const Vector<T>& y);

		// euclidean norm : norm2(A) =sqrt( ps(A,A))
		template<typename T, typename MP, typename MS, typename MO>
		T							norm2(const Matrix<T, MP, MS, MO>& x);

		template<typename T>
		T							norm2(const Vector<T>& x);

		// ritorna l'indice vettoriale della componente di valore massimo dell'Matrix (1-base)
		template<typename T, typename MP, typename MS, typename MO>
		size_t						imax(const Matrix<T, MP, MS, MO>& a);

		template<typename T>
		size_t						imax(const Vector<T>& a);

		// ritorna l'indice vettoriale della componente di valore minimo dell'Matrix (1-base)
		template<typename T, typename MP, typename MS, typename MO>
		size_t						imin(const Matrix<T, MP, MS, MO>& a);
		template<typename T>
		size_t						imin(const Vector<T>& a);



		// matrix product GE * GE (BLAS GEMM)
		// C= op(A) * op(B)      where  op(A)=A if  trA==false and  op(A)=A^T if trA==true
		template<typename T, typename MP, typename MS, typename MO>
		Matrix<T, MP, MS, MO>&		product(const	Matrix<T, MP, MS, MO>& A,
									 		const	Matrix<T, MP, MS, MO>& B,
													Matrix<T, MP, MS, MO>& C,
											bool trA = false, bool trB = false);


		template<typename T, typename MP, typename MS, typename MO>
		Matrix<T, MP, MS, MO>&		product(const	Matrix<T, MP, MS, MO>& A,
											const	Vector<T>& B,
													Vector<T>& C,
												bool trA = false);



		/// matrix product  GE * SYM (BLAS SYMM)
		template<typename T, typename MO>
		Matrix<T, M_PROP::GE, M_SHAPE::RE, MO>& product(const	Matrix<T, M_PROP::GE, M_SHAPE::RE, MO>& A,
														const	Matrix<T, M_PROP::SYM, M_SHAPE::RE, MO>& B,
																Matrix<T, M_PROP::GE, M_SHAPE::RE, MO>& C);


		/// matrix product   SYM * GE (BLAS SYMM)
		template<typename T, typename MO>
		Matrix<T, M_PROP::GE, M_SHAPE::RE, MO>& product(const	Matrix<T, M_PROP::SYM, M_SHAPE::RE, MO>& A,
														const	Matrix<T, M_PROP::GE, M_SHAPE::RE, MO>& B,
																Matrix<T, M_PROP::GE, M_SHAPE::RE, MO>& C);


		template<typename T, typename MO>
		Matrix<T, M_PROP::GE, M_SHAPE::RE, MO>& product(const	Matrix<T, M_PROP::SYM, M_SHAPE::RE, MO>& A,
														const	Vector<T>& B,	Vector<T>& C);




		// C=diag{A1,A2,...,AN} * B
		// note that diagonal block matrices must be passed by pointers
		// e.g., if A1,A2,..., B,C are matrices then use
		// diag_block_mat_product({&A1,&A2,...},B,C);
		template<typename T, typename MP, typename MS, typename MO>
		Matrix<T, MP, MS, MO>&		diag_block_mat_product(std::initializer_list< const Matrix<T, MP, MS, MO> *> A, const Matrix<T, MP, MS, MO>& B, Matrix<T, MP, MS, MO>& C);







		////////////////////////////////////////////////////////////////////////////////////////////////////
		//										IMPLEMENTATION
		////////////////////////////////////////////////////////////////////////////////////////////////////


		//-------------------------------
		// y := y + a x
		//-------------------------------
		template<typename T, typename MP, typename MS, typename MO>
		inline Matrix<T, MP, MS, MO>& axpy(T a, const Matrix<T, MP, MS, MO>& x, Matrix<T, MP, MS, MO>& y)
		{
			if (x.nrows() != y.nrows() || x.ncols() != y.ncols())
				throw std::range_error("Matrix& axpy(T a, const Matrix& x,  Matrix& y): size mismatch");

			Cblas<T>::axpy((const int)x.size(), a, x.data(), y.data());
			return y;
		}

		template<typename T>
		inline Vector<T>& axpy(T a, const Vector<T>& x, Vector<T>& y)
		{
			if (x.size() != y.size())
				throw std::range_error("Vector& axpy: size mismatch");
			Cblas<T>::axpy((const int)x.size(), a, x.data(), y.data());
			return y;
		}

		//-------------------------------
		//  T:= x*y : scalar product between vectors
		//-------------------------------
		template<typename T, typename MP, typename MS, typename MO>
		inline T dot(const Matrix<T, MP, MS, MO>& x, const Matrix<T, MP, MS, MO>& y)
		{
			size_t N = size_t(x.size());
			if (N != size_t(y.size()))
				throw std::range_error("T MML3::scalar_product(const Matrix<T>& x, const Matrix<T>& y): size mismatch");

			return Cblas<T>::dot((int)N, x.data(), y.data());
		}
		template<typename T>
		inline T dot(const Vector<T>& x, const Vector<T>& y)
		{

			if (x.size() != y.size())
				throw std::range_error(" vector dot : size mismatch");

			return Cblas<T>::dot(x.size(), x.data(), 1,  y.data(), 1);
		}

		//-------------------------------
		//  T:= || x || : euclidean norm
		//-------------------------------
		template<typename T, typename MP, typename MS, typename MO>
		inline T  norm2(const Matrix<T, MP, MS, MO>& x)
		{
			return Cblas<T>::nrm2((const int)x.size(), x.data());
		}

		template<typename T>

		inline T  norm2(const Vector<T>& x)
		{
			return Cblas<T>::nrm2(x.size(), x.data());
		}


		//-------------------------------
		//  x:= val * x : vector scaling
		//-------------------------------
		template<typename T, typename MP, typename MS, typename MO>
		inline Matrix<T, MP, MS, MO>& scal(Matrix<T, MP, MS, MO>& x, T val)
		{
			Cblas<T>::scal(x.size(), val, x.data());
			return x;
		}

		template<typename T>
		inline Vector<T>& scal(Vector<T>& x, T val)
		{
			Cblas<T>::scal(x.size(), val, x.data());
			return x;
		}


		//-------------------------------
		//  i:= arg_max(x(i))
		//-------------------------------
		template<typename T, typename MP, typename MS, typename MO>
		inline size_t imax(const Matrix<T, MP, MS, MO>& a)
		{
			return Cblas<T>::iamax((const int)a.size(), a.data(), 1) + MML3::BASE::OFFSET;
		}

		template<typename T>
		inline size_t imax(const Vector<T>& a)
		{
			return Cblas<T>::iamax((const int)a.size(), a.data(), 1) + MML3::BASE::OFFSET;
		}
		//-------------------------------
		//  i:= arg_min(x(i))
		//-------------------------------
		template<typename T, typename MP, typename MS, typename MO>
		inline size_t imin(const Matrix<T, MP, MS, MO>& a)
		{
			return Cblas<T>::iamin((const int)a.size(), a.data(), 1) + MML3::BASE::OFFSET;
		}

		template<typename T>
		inline size_t imin(const Vector<T>& a)
		{
			return Cblas<T>::iamin((const int)a.size(), a.data(), 1) + MML3::BASE::OFFSET;
		}



		//---------------------------------------------
		//  C= op(A) * op(B)  with A,B General matrices
		//---------------------------------------------
		template<typename T, typename MO>
		Matrix<T, M_PROP::GE, M_SHAPE::RE, MO>& product(const	Matrix<T, M_PROP::GE, M_SHAPE::RE, MO>& A,
														const	Matrix<T, M_PROP::GE, M_SHAPE::RE, MO>& B,
																Matrix<T, M_PROP::GE, M_SHAPE::RE, MO>& C=B,
														bool trA = false, bool trB = false)
		{
			bool tmp = false;
			

			
			Option::Transpose TRA = trA ? Option::Trans : Option::NoTrans;
			Option::Transpose TRB = trB ? Option::Trans : Option::NoTrans;

			// if needed resize C
			size_t M = trA ? A.ncols() : A.nrows();
			size_t N = trB ? B.nrows() : B.ncols();
			
			size_t K = trA ? A.nrows() : A.ncols();
			size_t K1 = trB ? B.ncols() : B.nrows();
			if (K1 != K)
				throw std::runtime_error("Cbla::gemm size");
			
			
			if (&B == &C)
			{
				Matrix<T, M_PROP::GE, M_SHAPE::RE, MO> Btmp(B);
				C.resize(M, N);
				Cblas<T>::gemm(A.CblasOrder(), TRA, TRB, M, N, K, T(1), A.data(), A.leading_dim(), Btmp.data(), Btmp.leading_dim(), T(0), C.data(), C.leading_dim());
			}
			else
			{
				C.resize(M, N);
				Cblas<T>::gemm(A.CblasOrder(), TRA, TRB, M, N, K, T(1), A.data(), A.leading_dim(), B.data(), B.leading_dim(), T(0), C.data(), C.leading_dim());
			}

			return C;
		}

		template<typename T, typename MO>
		Vector<T>& product(	const	Matrix<T, M_PROP::GE, M_SHAPE::RE, MO>& A,
							const	Vector<T>& B,
									Vector<T>& C,
									bool trA = false)
		{
			bool tmp = false;
			const Vector<T>* pB = &B;

			if (&B == &C)
			{
				pB = new  Vector<T>(B);
				tmp = true;
			}
			Option::Transpose TRA = trA ? Option::Trans : Option::NoTrans;


			// if needed resize C
			size_t szC = trA ? A.ncols() : A.nrows();
			C.resize(szC);
			size_t K = trA ? A.nrows() : A.ncols();

			if (B.size()!= K)
				throw std::runtime_error("Cbla::gemm size");
			Cblas<T>::gemv(A.CblasOrder(), TRA, A.nrows(), A.ncols(), T(1), A.data(), A.leading_dim(), pB->data(), 1, T(0), C.data(), 1);
			if (tmp)
				delete pB;
			return C;
		}
		

			//---------------------------------------------
			//  C= A * B  with A General and B symmetric
			//---------------------------------------------
			template<typename T, typename MO>
			Matrix<T, M_PROP::GE, M_SHAPE::RE, MO>& product(const	Matrix<T, M_PROP::GE,	M_SHAPE::RE, MO>& A,
															const	Matrix<T, M_PROP::SYM,	M_SHAPE::RE, MO>& B,
																	Matrix<T, M_PROP::GE,	M_SHAPE::RE, MO>& C)
			{
				// test di congruità dimensionale
				if (A.ncols() != B.nrows())
					throw std::runtime_error("MML3::Algorithm:: product( Matrix , symMatrix,...): incompatible matrix dimensions");

				index_type M = A.nrows(), N = B.ncols();
				if (C.nrows() != M || C.ncols() != N)
					C.resize(M, N);

				Cblas<T>::symm(B.CblasOrder(), Option::Right,
					B.CblasSymUpLo(), M, N,
					T(1), B.data(), B.leading_dim(),
					A.data(), A.leading_dim(), T(0),
					C.data(), C.leading_dim());
				return C;
			}


			//---------------------------------------------
			//  C= A * B  with A Symmetric and B GEneral
			//---------------------------------------------
			template<typename T, typename MO>
			Matrix<T, M_PROP::GE, M_SHAPE::RE, MO>& product(const Matrix<T, M_PROP::SYM, M_SHAPE::RE, MO>& A,
															const Matrix<T, M_PROP::GE,  M_SHAPE::RE, MO>&  B,
																  Matrix<T, M_PROP::GE,  M_SHAPE::RE, MO>& C)
			{
				// test di congruità dimensionale
				if (A.ncols() != B.nrows())
					throw std::runtime_error("MML3::ALG:: product( symMatrix , Matrix,...): incompatible matrix dimensions");
				index_type M = A.nrows(), N = B.ncols();
				if (C.nrows() != M || C.ncols() != N)
					C.resize(M, N);
				Cblas<T>::symm(A.CblasOrder(), Option::Left,
					A.CblasSymUpLo(), M, N,
					T(1), A.data(), A.leading_dim(),
					B.data(), B.leading_dim(), T(0),
					C.data(), C.leading_dim());
				return C;
			}








		/*
		template<typename T>
		Matrix<T, typename M_PROP::GE>& diag_block_mat_product(std::initializer_list< const Matrix<T, M_PROP::GE> *> A, const Matrix<T, M_PROP::GE>& B, Matrix<T, M_PROP::GE>& C)
		{

			// test di conformita' delle dimensioni
			size_t szA = 0;
			for (auto a : A)
			{
				if (a->nrows() != a->ncols())
					throw std::length_error("MML3::diag_block_mat_product: all diagonal blocks must be square matrices");
				szA += a->nrows();
			}
			if (szA != B.nrows())
				throw std::length_error("MML3::diag_block_mat_product: size of A is incompatible with the size of B");

			// prodotto
			C.resize(szA, B.ncols()).fill(T(0));

			// indica il primo indice diagonale della sottomatrice in A ( meno 1)
			size_t sz = 0;
			size_t Bcols = B.ncols();

			// ciclo sui blocchi diagonali di A
			for (auto a : A)
			{
				const Matrix<T>& aa = *a; // alias per il blocco diagonale
				size_t sz_aa = aa.ncols();
				size_t I0 = sz;
				for (size_t i = 0; i != sz_aa; ++i)
				{
					T*			C_I = C[I0 + i];
					const T*	a_i = aa[i];
					// cycles swapped for efficiency
					for (size_t k = 0; k != sz_aa; ++k)
					{
						T a_ik = a_i[k];
						const T* B_K = B[I0 + k];
						for (size_t J = 0; J != Bcols; ++J)
							C_I[J] += a_ik * B_K[J];
					}
				}
				sz += sz_aa;
			}
			return C;

		}


		*/





	}// end namespace ALG
} // end namespace MML3
