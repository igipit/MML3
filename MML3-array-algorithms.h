#pragma once
#include<algorithm>
#include<type_traits>

namespace MML3{

#define FORALL(i,n) for (SZ i = 0; i != n; ++i)

	template<typename T>
	struct ArrayOp
	{

		static_assert(std::is_arithmetic<T>::value,	"MML3::ArrayOp requires arithmetic elements");

		typedef size_t		SZ;
		typedef T*			P;
		typedef const T*	const_P;

		



		static P copy(SZ sz, const_P x, P z)					{ memcpy(z, x, sz*sizeof(T)); return z; }

		static P fill(SZ sz, T val, P z)						{ std::fill_n(z, sz, val); return z; }

		static P zero(SZ sz, P z)								{ memset(z, 0, sz*sizeof(T)); return z; }

		static P scal(SZ sz, T scaler, P z)						{ FORALL(i, sz)	z[i] *= scaler; return z; }

		static P sum(SZ sz, const_P x, const_P y, P z)			{ FORALL(i, sz)	z[i] = x[i] + y[i]; return z; }

		static P sub(SZ sz, const_P x, const_P y, P z)			{ FORALL(i, sz)	z[i] = x[i] - y[i]; return z; }

		static P axpy(SZ sz, T a, const_P x, P y)
		{ 
			if (a == T(1)) 
				FORALL(i, sz)	y[i] += x[i];
			else if (a == T(-1))
				FORALL(i, sz)	y[i] -= x[i];
			else
				FORALL(i, sz)	y[i] += a*x[i]; 
			return y;
		}



		static T	dot(SZ sz, const_P x, const_P y)				{ T acc = T(0);   FORALL(i, sz)	acc += x[i] * y[i]; return acc; }

		static T	norm(SZ sz, const_P x)							{ return sqrt(dot(sz, x, x)); }

		static T	unitize(SZ sz, P z)								{ T tmp = norm(sz, z); scal(sz, 1 / tmp, z); return tmp; }

		
		static P lin_comb(SZ sz, T a1, const_P x1, T a2, const_P x2, P z)
		{
			FORALL(i, sz)
				z[i] = a1*x1[i] + a2* x2[i];
			return z;
		}

		static P lin_comb(SZ sz, T a1, const_P x1, T a2, const_P x2, T a3, const_P x3, P z)
		{
			FORALL(i, sz)
				z[i] = a1*x1[i] + a2* x2[i] + a3* x3[i];
			return z;
		}

		static P lin_comb(SZ sz, SZ n, const T a[], const_P x[], P z)
		{
			zero(sz,  z);
			FORALL(i, n) 	axpy(sz, a[i], x[i], z);
			return z;
		}
	};

#undef FORALL




	template<int N, typename T = double>
	struct StaticArrayOp
	{

		

		typedef T*			P;
		typedef const T*	const_P;

		static_assert(std::is_arithmetic<T>::value, "MML3::StaticArrayOp requires arithmetic elements");


		static  P fill(T value, P z){return ArrayOp<double>::fill(N, value, z);	}

		static  P zero(P z)	{ return ArrayOp<double>::zero(N, z);
		}


		static  P scal(T value, P z){return ArrayOp<double>::scal(N, value, z);	}

		// prodotto scalare tra vettori
		static  T dot(const_P a, const_P b)	{return ArrayOp<double>::dot(N, a, b);}
		// norma euclidea di un vettore
		static inline T norm2(const_P a){return ArrayOp<double>::norm(N, a);}




		static inline T unitize(P z){return ArrayOp<double>::unitize(N, z);	}

		static inline P lin_comb(T a1, const_P x1, T a2, const_P x2, P z){ return 	ArrayOp<double>::lin_comb(N, a1, x1, a2, x2, z);}

		static inline P lin_comb(int n, const T a[], const_P x[], P z) { return ArrayOp<double>::lin_comb(N, n, a, x, z);}

		static inline P axpy(const T a, const_P x, P y)	{	return ArrayOp<double>::axpy(N, a, x, y);}

		/*static inline void add(const_P x, P y)
		{
			for (auto i = 0; i != N; ++i)
				y[i] += x[i];
		}*/

		static inline P sum(const_P x, const_P y, P z){return 	ArrayOp<double>::sum(N, x, y, z);}

		static inline P sub(const_P x, const_P y, P z){ return ArrayOp<double>::sub(N, x, y, z);	}




		static inline P copy(const_P source, P dest){ return ArrayOp<double>::copy(N, source, dest);}
		

		// prodotto vettore 3D

		static inline P cross_product(const_P a, const_P b, P c)
		{

			static_assert(N == 3, "MML3::StaticArrayOp cross product defined only for 3D vectors");
			c[0] = a[1] * b[2] - a[2] * b[1];
			c[1] = a[2] * b[0] - a[0] * b[2];
			c[2] = a[0] * b[1] - a[1] * b[0];
			return c;
		}

		// prodotto vettore2D
		static inline T cross_product(const_P a, const_P b)
		{
			static_assert(N == 2, "MML3::StaticArrayOp  cross product defined only for 2D vectors");
			return a[0] * b[1] - a[1] * b[0];
		}


		
	};



}// end MML3