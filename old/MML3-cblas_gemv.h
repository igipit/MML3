#ifndef CBLAS_GEMV_H_INCLUDED
#define CBLAS_GEMV_H_INCLUDED
#include"MML3-config.h"
namespace AnOnYmOuS
{
template<typename I>
inline size_t OFFSET(I N, int incX){ return (incX > 0) ? 0 : (N-1)*(-incX); }
}


// y= beta*y + alpha A*x or y= beta*y + alpha A'*x
template<typename T>
void cblas_gemv(MML3::CblasOption::Order order,
                MML3::CblasOption::Transpose TransA,
                size_t M,   size_t N, T alpha, const T *A,
                size_t lda, T *X, MML3::lapack_int_t incX,
                T beta, T *Y, MML3::lapack_int_t incY)
{

    using AnOnYmOuS::OFFSET;
    static_assert(std::is_same<T,float>::value ||std::is_same<T,double>::value," only for reals" );
    size_t i, j;
    size_t lenX, lenY;

    bool Trans = !(TransA==MML3::CblasOption::Transpose::NoTrans);

    bool isRowMajor=(order==MML3::CblasOption::Order::RowMajor);

 //CHECK_ARGS12(GEMV,order,TransA,M,N,alpha,A,lda,X,incX,beta,Y,incY);

      if (M == 0 || N == 0)
        return;

      if (alpha == 0.0 && beta == 1.0)
        return;

      if (!Trans)
        {
        lenX = N;
        lenY = M;
      } else {
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
      } else if (beta != 1.0) {
        size_t iy = OFFSET(lenY, incY);
        for (i = 0; i < lenY; i++) {
          Y[iy] *= beta;
          iy += incY;
        }
      }

      if (alpha == 0.0)
        return;

      if ( (isRowMajor && !Trans) || (!isRowMajor && Trans) )
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
      else if ( (isRowMajor && Trans )  || (!isRowMajor && !Trans) )
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
      } else
      {
        throw std::runtime_error("cblas_gemv: unrecognized operation");
      }
    }


#endif // CBLAS_GEMV_H_INCLUDED
