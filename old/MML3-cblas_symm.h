

#include"MML3-config.h"

template<typename T>
void cblas_symm	(
	MML3::CblasOption::Order 	order_,
	MML3::CblasOption::Side 	side_,
	MML3::CblasOption::UpLo 	uplo_,
	size_t	M,size_t 	N,	T 	alpha,
	const T * 	A, size_t 	lda,
	const T * 	B, size_t 	ldb,
	 T 	beta,T * 	C,size_t 	ldc )
{
	size_t i, j, k;
	size_t n1, n2;

	bool isRowMajor=(order_ == MML3::CblasOption::Order::RowMajor);
	bool isUpper=(uplo_==MML3::CblasOption::UpLo::Upper);
	bool isLeft=(side_==MML3::CblasOption::Side::Left);


      //CHECK_ARGS13(SYMM,Order,Side,Uplo,M,N,alpha,A,lda,B,ldb,beta,C,ldc);

      if (alpha == 0.0 && beta == 1.0)
        return;

      if (isRowMajor) {
        n1 = M;
        n2 = N;
      } else {
        n1 = N;
        n2 = M;
        isUpper=!isUpper;
		isLeft=!isLeft;
      }

      /* form  y := beta*y */
      if (beta == 0.0) {
        for (i = 0; i < n1; i++) {
          for (j = 0; j < n2; j++) {
            C[ldc * i + j] = 0.0;
          }
        }
      } else if (beta != 1.0) {
        for (i = 0; i < n1; i++) {
          for (j = 0; j < n2; j++) {
            C[ldc * i + j] *= beta;
          }
        }
      }

      if (alpha == 0.0)
        return;

      if (isLeft && isUpper) {

        /* form  C := alpha*A*B + C */

        for (i = 0; i < n1; i++) {
          for (j = 0; j < n2; j++) {
            const T temp1 = alpha * B[ldb * i + j];
            T temp2 = 0.0;
            C[i * ldc + j] += temp1 * A[i * lda + i];
            for (k = i + 1; k < n1; k++) {
              const T Aik = A[i * lda + k];
              C[k * ldc + j] += Aik * temp1;
              temp2 += Aik * B[ldb * k + j];
            }
            C[i * ldc + j] += alpha * temp2;
          }
        }

      } else if (isLeft && !isUpper) {

        /* form  C := alpha*A*B + C */

        for (i = 0; i < n1; i++) {
          for (j = 0; j < n2; j++) {
            const T temp1 = alpha * B[ldb * i + j];
            T temp2 = 0.0;
            for (k = 0; k < i; k++) {
              const T Aik = A[i * lda + k];
              C[k * ldc + j] += Aik * temp1;
              temp2 += Aik * B[ldb * k + j];
            }
            C[i * ldc + j] += temp1 * A[i * lda + i] + alpha * temp2;
          }
        }

      } else if (!isLeft && isUpper) {

        /* form  C := alpha*B*A + C */

        for (i = 0; i < n1; i++) {
         for (j = 0; j < n2; j++) {
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

     } else if (!isLeft && !isUpper) {

       /* form  C := alpha*B*A + C */

       for (i = 0; i < n1; i++) {
         for (j = 0; j < n2; j++) {
           const T temp1 = alpha * B[ldb * i + j];
           T temp2 = 0.0;
           for (k = 0; k < j; k++) {
             const T Ajk = A[j * lda + k];
             C[i * ldc + k] += temp1 * Ajk;
             temp2 += B[ldb * i + k] * Ajk;
           }
           C[i * ldc + j] += temp1 * A[j * lda + j] + alpha * temp2;
         }
       }

     } else
	 {
       throw std::runtime_error("cblas_symm: unrecognized operation");
     }
   }
