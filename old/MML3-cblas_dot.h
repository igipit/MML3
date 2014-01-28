
#pragma once
namespace MML3
{
	namespace cblas
	{
		template<typename T>
		T	cblas_sdot(size_t N, const T* X, int incX, const T* Y, int incY)
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
	} // end namespace
} // end namespace