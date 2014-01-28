
   #include"MML3-config.h"

namespace MML3
{
	namespace my_cblas
	{

		template<typename T>
		void cblas_gemm(CblasOption::Order Order,
			CblasOption::Transpose Trans_A,
			CblasOption::Transpose Trans_B,
			size_t M, size_t N, size_t K,
			T alpha,
			const T *A, size_t lda,
			const T *B, size_t ldb,
			T beta,
			T	 *C, size_t ldc)
		{

			static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value, " only for reals");
			size_t i, j, k;
			size_t n1, n2;
			size_t ldf, ldg;
			bool TransF, TransG;

			bool TransA = !(Trans_A == CblasOption::Transpose::NoTrans);
			bool TransB = !(Trans_B == CblasOption::Transpose::NoTrans);
			const T *F, *G;

			//CHECK_ARGS14(GEMM,Order,TransA,TransB,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc);

			if (alpha == 0.0 && beta == 1.0)
				return;

			if (Order == CblasOption::RowMajor)
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
				for (i = 0; i < n1; i++)
				for (j = 0; j < n2; j++)
					C[ldc * i + j] *= beta;
			}

			if (alpha == 0.0)
				return;

			if (!TransF && !TransG)
			{
				/* form  C := alpha*A*B + C */
				for (k = 0; k < K; k++)
				for (i = 0; i < n1; i++) {
					const T temp = alpha * F[ldf * i + k];
					if (temp != 0.0)
					{
						for (j = 0; j < n2; j++)
							C[ldc * i + j] += temp * G[ldg * k + j];
					}
				}

			}
			else if (!TransF &&   TransG)
			{
				/* form  C := alpha*A*B' + C */
				for (i = 0; i < n1; i++)
				for (j = 0; j < n2; j++)
				{
					T temp = 0.0;
					for (k = 0; k < K; k++)
						temp += F[ldf * i + k] * G[ldg * j + k];
					C[ldc * i + j] += alpha * temp;
				}
			}
			else if (TransF &&   !TransG)
			{
				for (k = 0; k < K; k++)
				for (i = 0; i < n1; i++) {
					const T temp = alpha * F[ldf * k + i];
					if (temp != 0.0) {
						for (j = 0; j < n2; j++) {
							C[ldc * i + j] += temp * G[ldg * k + j];
						}
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
				throw std::runtime_error("cblas_gemm: unrecognized operation");
			}
		}


	}// end namespace
}// end namespace
