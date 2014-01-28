#pragma once
#include"MML3-config.h"


// ===========================================================================
//
// C++ cblas implementation for some cblas functions
// ===========================================================================
namespace MML3{





	template< typename T>
	struct Cblas
	{

		// functions implemented only for reals
		static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value, " only for reals");

		///@name  Cblas LEVEL 1
		//@{

		//  sum_i X_i * Y_i
		static T		dot(const size_t N, const T  *X, const int incX, const T  *Y, const int incY);
		static   T		nrm2(const size_t N, const T* X, const int incX = 1);

		// Y <--alpha X + Y
		static void		axpy(const size_t N, const T alpha, const T *X, const int incX, T *Y, const int incY);
			// X <-- alpha*X
		static void		scal(const size_t N, const T alpha, T *X, const int incX = 1);
		static size_t	iamax(const size_t N, const T  *X, const int incX);
		static size_t	iamin(const size_t N, const T  *X, const int incX);
		//@}

		///@name  Cblas LEVEL 2
		//@{

		static void gemv(const Option::Order orderA, const  Option::Transpose TransA, const size_t M, const size_t N, const T alpha, const T *A, const size_t lda, const T *X, const int incX, const T beta, T *Y, const int incY);



		//@}

		///@name  Cblas LEVEL 3
		//@{
		static void gemm(const  Option::Order Order, const  Option::Transpose TransA, const  Option::Transpose TransB,
			const size_t M, const size_t N, const size_t K, const T alpha,	const T *A, const int lda,	const T *B, const int ldb,
			const T beta,T	 *C, const int ldc);
		//@}



		static void symm(const  Option::Order Order, const  Option::Side Side,
			const  Option::UpLo Uplo, const size_t M, const size_t N,
			const T alpha, const T *A, const int lda,const T *B, const int ldb, const T beta,
			T *C, const int ldc);

	};





	template<typename T>
	inline T	Cblas<T>::dot(const size_t N, const T* X, const int incX, const T* Y, const int incY)
	{
		T tmp = T(0);
		if (incX == 1 && incY == 1)
		{
			for (size_t i = 0; i != N; ++i)
				tmp += X[i] * Y[i];
		}
		else
		{
			for (size_t i = 0; i != N; ++i)
			{
				tmp += (*X) * (*Y);
				X += incX;
				Y += incY;
			}
		}
		return tmp;
	}

	template<typename T>
	inline T	Cblas<T>::nrm2(const size_t N, const T* X, const int incX)
	{
		return sqrt(dot(N, X, incX, X, incX));
	}


	template<typename T>
	inline void Cblas<T>::axpy(const size_t N, const T alpha, const T *X, const int incX, T *Y, const int incY)
	{
		if (alpha == T(0))
			return;
		if (incX == 1 && incY == 1)
		{
			for (size_t i = 0; i != N; ++i)
				Y[i] += alpha*X[i];
		}
		else
		{
			for (size_t i = 0; i != N; ++i)
			{
				*Y += alpha * (*X);
				X += incX;
				Y += incY;
			}
		}
	}

	template<typename T>
	inline void Cblas<T>::scal(const size_t N, const T alpha, T *X, const int incX)
	{
		if (alpha == T(1))
			return;
		if (incX == 1)
		{
			for (size_t i = 0; i != N; ++i)
				X[i] = alpha*X[i];
		}
		else
		{
			for (size_t i = 0; i != N; ++i)
			{
				*X += alpha * (*X);
				X += incX;
			}
		}
	}


	template<typename T>
	inline size_t Cblas<T>::iamax(const size_t N, const T  *X, const int incX)
	{
		if (!N)
			return 0;
		T		max_val = X[0];
		size_t	max_idx = 0;
		X += incX;
		for (size_t i = 1; i != N; ++i)
		{
			if (*X > max_val)
			{
				max_idx = i;
				max_val = *X;
			}
			X += incX;
		}
		return max_idx;
	}

	template<typename T>
	inline size_t Cblas<T>::iamin(const size_t N, const T  *X, const int incX)
	{
		if (!N)
			return 0;
		T		min_val = X[0];
		size_t	min_idx = 0;
		X += incX;
		for (size_t i = 1; i != N; ++i)
		{
			if (*X < min_val)
			{
				min_idx = i;
				min_val = *X;
			}
			X += incX;
		}
		return min_idx;
	}


	namespace AnOnYmOuS
	{
		template<typename I>
		inline size_t OFFSET(I N, int incX){ return (incX > 0) ? 0 : (N - 1)*(-incX); }
	}



	template<typename T>
	inline void Cblas<T>::gemv(const Option::Order orderA, const  Option::Transpose TransA, const size_t M, const size_t N,
		const T alpha, const T *A, const size_t lda, const T *X, const int incX, const T beta, T *Y, const int incY)
	{

		using AnOnYmOuS::OFFSET;

		size_t i, j;
		size_t lenX, lenY;

		bool Trans = !(TransA == Option::Transpose::NoTrans);

		bool isRowMajor = (orderA == Option::Order::RowMajor);

		//CHECK_ARGS12(GEMV,order,TransA,M,N,alpha,A,lda,X,incX,beta,Y,incY);

		if (M == 0 || N == 0)
			return;

		if (alpha == 0.0 && beta == 1.0)
			return;

		if (!Trans)
		{
			lenX = N;
			lenY = M;
		}
		else {
			lenX = M;
			lenY = N;
		}

		/* form  y := beta*y */
		if (beta == 0.0) {
			size_t iy = OFFSET(lenY, incY);
			for (i = 0; i < lenY; i++) {
				Y[iy] = 0.0;
				iy += incY;
			}
		}
		else if (beta != 1.0) {
			size_t iy = OFFSET(lenY, incY);
			for (i = 0; i < lenY; i++) {
				Y[iy] *= beta;
				iy += incY;
			}
		}

		if (alpha == 0.0)
			return;

		if ((isRowMajor && !Trans) || (!isRowMajor && Trans))
		{
			/* form  y := alpha*A*x + y */
			size_t iy = OFFSET(lenY, incY);
			for (i = 0; i < lenY; i++)
			{
				T temp = 0.0;
				size_t ix = OFFSET(lenX, incX);
				for (j = 0; j < lenX; j++) {
					temp += X[ix] * A[lda * i + j];
					ix += incX;
				}
				Y[iy] += alpha * temp;
				iy += incY;
			}
		}
		else if ((isRowMajor && Trans) || (!isRowMajor && !Trans))
		{
			/* form  y := alpha*A'*x + y */
			size_t ix = OFFSET(lenX, incX);
			for (j = 0; j < lenX; j++)
			{
				const T temp = alpha * X[ix];
				if (temp != 0.0)
				{
					size_t iy = OFFSET(lenY, incY);
					for (i = 0; i < lenY; i++)
					{
						Y[iy] += temp * A[lda * j + i];
						iy += incY;
					}
				}
				ix += incX;
			}
		}
		else
		{
			assert(false);
		}
	}



	template<typename T>
	inline void Cblas<T>::gemm(const  Option::Order Order, const  Option::Transpose Trans_A, const  Option::Transpose Trans_B,
		const size_t M, const size_t N, const size_t K,
		const T alpha,
		const T *A, const int lda,
		const T *B, const int ldb,
		const T beta,
		T	 *C, const int ldc)
	{

		size_t i, j, k;
		size_t n1, n2;
		size_t ldf, ldg;
		bool TransF, TransG;

		bool TransA = !(Trans_A == Option::Transpose::NoTrans);
		bool TransB = !(Trans_B == Option::Transpose::NoTrans);
		const T *F, *G;

		if (alpha == 0.0 && beta == 1.0)
			return;

		if (Order == Option::RowMajor)
		{
			n1 = M;
			n2 = N;
			F = A;
			ldf = lda;
			TransF = TransA;
			G = B;
			ldg = ldb;
			TransG = TransB;
		}
		else
		{
			n1 = N;
			n2 = M;
			F = B;
			ldf = ldb;
			TransF = TransB;
			G = A;
			ldg = lda;
			TransG = TransA;
		}

		/* form  y := beta*y */
		if (beta == 0.0) {
			for (i = 0; i != n1; i++)
			for (j = 0; j != n2; j++)
				C[ldc * i + j] = 0.0;

		}
		else if (beta != 1.0)
		{
			for (i = 0; i != n1; i++)
			for (j = 0; j != n2; j++)
				C[ldc * i + j] *= beta;
		}

		if (alpha == 0.0)
			return;


		if (!TransF && !TransG)
		{
			/* form  C := alpha*A*B + C */
			for (k = 0; k != K; k++)
			{
				for (i = 0; i != n1; i++)
				{
					const T temp = alpha * F[ldf * i + k];
					if (temp != 0.0)
					{
						for (j = 0; j != n2; j++)
							C[ldc * i + j] += temp * G[ldg * k + j];

					}
				}
			}

		}
		else if (!TransF &&   TransG)
		{
			/* form  C := alpha*A*B' + C */
			for (i = 0; i !=n1; i++)
			for (j = 0; j !=n2; j++)
			{
				T temp = 0.0;
				for (k = 0; k != K; k++)
					temp += F[ldf * i + k] * (G[ldg * j + k]);
				C[ldc * i + j] += alpha * temp;
			}
		}
		else if (TransF &&   !TransG)
		{
			for (k = 0; k < K; k++)
			for (i = 0; i < n1; i++)
			{
				const T temp = alpha * F[ldf * k + i];
				if (temp != 0.0)
				{
					for (j = 0; j < n2; j++)
						C[ldc * i + j] += temp * G[ldg * k + j];

				}
			}

		}
		else if (TransF && TransG)
		{
			for (i = 0; i < n1; i++)
			for (j = 0; j < n2; j++)
			{
				T temp = 0.0;
				for (k = 0; k < K; k++)
					temp += F[ldf * k + i] * G[ldg * j + k];
				C[ldc * i + j] += alpha * temp;
			}
		}
		else
		{
			assert(false);
		}


	}


	template<typename T>
	inline void Cblas<T>::symm(const  Option::Order Order, const  Option::Side Side,
		const  Option::UpLo Uplo, const size_t M, const size_t N,
		const T alpha, const T *A, const int lda,
		const T *B, const int ldb, const T beta,
		T *C, const int ldc)
	{


		size_t i, j, k;
		size_t n1, n2;

		bool isRowMajor = (Order == MML3::Option::Order::RowMajor);
		bool isUpper = (Uplo == MML3::Option::UpLo::Upper);
		bool isLeft = (Side == MML3::Option::Side::Left);


		//CHECK_ARGS13(SYMM,Order,Side,Uplo,M,N,alpha,A,lda,B,ldb,beta,C,ldc);

		if (alpha == 0.0 && beta == 1.0)
			return;

		if (isRowMajor)
		{
			n1 = M;
			n2 = N;
		}
		else
		{
			n1 = N;
			n2 = M;
			isUpper = !isUpper;
			isLeft = !isLeft;
		}

		/* form  y := beta*y */
		if (beta == 0.0)
		{
			for (i = 0; i < n1; i++)
			{
				for (j = 0; j < n2; j++)
				{
					C[ldc * i + j] = 0.0;
				}
			}
		}
		else if (beta != 1.0)
		{
			for (i = 0; i < n1; i++)
			{
				for (j = 0; j < n2; j++)
				{
					C[ldc * i + j] *= beta;
				}
			}
		}

		if (alpha == 0.0)
			return;

		if (isLeft && isUpper)
		{
			/* form  C := alpha*A*B + C */
			for (i = 0; i < n1; i++)
			{
				for (j = 0; j < n2; j++)
				{
					const T temp1 = alpha * B[ldb * i + j];
					T temp2 = 0.0;
					C[i * ldc + j] += temp1 * A[i * lda + i];
					for (k = i + 1; k < n1; k++)
					{
						const T Aik = A[i * lda + k];
						C[k * ldc + j] += Aik * temp1;
						temp2 += Aik * B[ldb * k + j];
					}
					C[i * ldc + j] += alpha * temp2;
				}
			}

		}
		else if (isLeft && !isUpper)
		{
			/* form  C := alpha*A*B + C */
			for (i = 0; i < n1; i++)
			{
				for (j = 0; j < n2; j++)
				{
					const T temp1 = alpha * B[ldb * i + j];
					T temp2 = 0.0;
					for (k = 0; k < i; k++)
					{
						const T Aik = A[i * lda + k];
						C[k * ldc + j] += Aik * temp1;
						temp2 += Aik * B[ldb * k + j];
					}
					C[i * ldc + j] += temp1 * A[i * lda + i] + alpha * temp2;
				}
			}

		}
		else if (!isLeft && isUpper)
		{
			/* form  C := alpha*B*A + C */
			for (i = 0; i < n1; i++)
			{
				for (j = 0; j < n2; j++)
				{
					const T temp1 = alpha * B[ldb * i + j];
					T temp2 = 0.0;
					C[i * ldc + j] += temp1 * A[j * lda + j];
					for (k = j + 1; k < n2; k++) {
						const T Ajk = A[j * lda + k];
						C[i * ldc + k] += temp1 * Ajk;
						temp2 += B[ldb * i + k] * Ajk;
					}
					C[i * ldc + j] += alpha * temp2;
				}
			}

		}
		else if (!isLeft && !isUpper)
		{

			/* form  C := alpha*B*A + C */

			for (i = 0; i < n1; i++)
			{
				for (j = 0; j < n2; j++)
				{
					const T temp1 = alpha * B[ldb * i + j];
					T temp2 = 0.0;
					for (k = 0; k < j; k++)
					{
						const T Ajk = A[j * lda + k];
						C[i * ldc + k] += temp1 * Ajk;
						temp2 += B[ldb * i + k] * Ajk;
					}
					C[i * ldc + j] += temp1 * A[j * lda + j] + alpha * temp2;
				}
			}

		}
		else
		{
			throw std::runtime_error("cblas_symm: unrecognized operation");
		}
	}




}  // end namespace MML3
