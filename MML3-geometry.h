#pragma once
#include<algorithm>
#include<vector>
#include"MML3-math.h"
#include"MML3-array-algorithms.h"


namespace MML3
{

	namespace GEO{



		// fixed size array to serve Point and Vector classes
		template<size_t N, typename CT = double>
		class array_
		{
			
		public:
			typedef  CT				REAL;
			typedef  CT*			Ptr;
			typedef const Ptr	const_Ptr;

			typedef MML3::StaticArrayOp<N, CT> VecOp;

			array_() = default;
			array_(const array_&) = default;
			array_(const std::initializer_list<REAL>& list){ std::copy_n(list.begin(), std::min(list.size(), N), x); }
			array_(const REAL v[]){ std::copy_n(v, N, x); }

			array_(REAL x1) :x{ x1 }{}
			array_(REAL x1, REAL x2) :x{ x1, x2 }{ static_assert(N > 1," bad dimension"); }
			array_(REAL x1, REAL x2,REAL x3) :x{ x1, x2,x3}{ static_assert(N > 2, " bad dimension"); }

			void copy(const REAL* src){ VecOp::copy(src, x); }

			explicit operator CT*(){ return x; }
			explicit operator const CT*()const{ return x; }
		
			array_& operator=(const array_&) = default;
			array_& operator=(const std::initializer_list<REAL>& list) { std::copy_n(list.begin(), std::min(list.size(), N), x); return *this; }



			REAL& operator()(size_t i){ return x[i - 1]; }
			const REAL& operator()(size_t i)const{ return x[i - 1]; }

			REAL& operator[](size_t i){ return x[i]; }
			const REAL& operator[](size_t i)const{ return x[i]; }

			size_t  size()const{ return N; }
			REAL*	data(){ return x; }
			const REAL*	data()const{ return x; }



			// algorithms

			REAL norm()const
			{
				return VecOp::norm2(x);	
			}

			array_& operator+=(const array_& o)
			{
				VecOp::sum(o.x, x, x);
				return *this;
			}

			array_& operator-=(const array_& o)
			{
				VecOp::sub(x,o.x,x);
				return *this;
			}

			array_& operator*=(REAL value)
			{
				VecOp::scal(value, x);
				return *this;
			}

			array_& operator/=(REAL value)
			{ 
				VecOp::scal(1/value, x); 
				return *this; 
			}

			REAL dot(const array_& o)const
			{
				return VecOp::dot(x, o.x);
			}

			array_& unitize()
			{ 
				VecOp::unitize(x);
				return *this; 
			}

		private:
			REAL x[N];

		}; // end class array_


		


		template<size_t N, typename CT = double>
		class Point : public  array_ < N, CT >
		{
			typedef  array_< N, CT > Array_;
		public:
	
			Point() = default;
			Point(const Point&) = default;
			Point(const Array_& o) :Array_(o){}

			Point(const std::initializer_list<REAL>& list) :Array_(list){}
			Point(const REAL v[]) :Array_(v){ }


			/*operator Array_&() { return *this; }
			operator const Array_&() const { return *this; }*/


			Point& operator=(const Point&) = default;
			Point& operator=(const Array_& o){ Array_::operator=(o); return *this; }
			Point& operator=(const std::initializer_list<REAL>& list){ Array_::operator=(list); return *this; }




			
		};



		template<size_t N, typename CT = double>
		class Vector : public  array_ < N, CT >
		{
			typedef  array_ < N, CT > Array_;
		public:

			Vector() = default;
			Vector(const Vector&) = default;
			Vector(const Array_& o) :Array_(o){}
			Vector(const std::initializer_list<REAL>& list) :Array_(list){}
			Vector(const REAL v[]) :Array_(v){ }

			Vector& operator=(const Vector&) = default;
			Vector& operator=(const Array_& o){ Array_::operator=(o); return *this; }
			Vector& operator=(const std::initializer_list<REAL>& list){ Array_::operator=(list); return *this; }
		};






		

		template<typename CT>
		array_<3, CT>& set_spherical(array_<3, CT>& P, CT r, CT latitude, CT longitude)
		{
			P[2] = r*cos(latitude);

			P[0] = P[2] * cos(longitude);
			P[1] = P[2] * sin(longitude);
			P[2] = r*sin(latitude);
			return P;
	
		}






		
		template<size_t N, typename CT>
		array_<N, CT>& axpy(CT a, const array_<N, CT>& x, array_<N, CT>& y)
		{
			array_<N, CT>::VecOp::axpy(a, x.data(), y.data());
			return y;
		}

		template<size_t N, typename CT>
		array_<N, CT>& lin_comb(CT a, const array_<N, CT>& x, CT b, const array_<N, CT>& y, array_<N, CT>& z)
		{

			array_<N, CT>::VecOp::lin_comb(a, x.data(), b, y.data(), z.data());
			return z;
		}

		template<size_t N, typename CT>
		array_<N, CT>& add(const array_<N, CT>& a, const array_<N, CT>& b, array_<N, CT>& dest)
		{
			array_<N, CT>::VecOp::sum(a.data(), b.data(), dest.data());
			return dest;
		}

		template<size_t N, typename CT>
		array_<N, CT>& sub(const array_<N, CT>& a, const array_<N, CT>& b, array_<N, CT>& dest)
		{
			array_<N, CT>::VecOp::sub(a.data(), b.data(), dest.data());
			return dest;
		}

		template< typename CT>
		array_<3, CT>  cross_product(const array_<3, CT>& x, const array_<3, CT>& y, array_<3, CT>& z)
		{
			CT *Z = z.data();
			const CT * a = x.data(), *b = y.data();

			Z[0] = a[1] * b[2] - a[2] * b[1];
			Z[1] = a[2] * b[0] - a[0] * b[2];
			Z[2] = a[0] * b[1] - a[1] * b[0];
			return z;

		}


		template<size_t N, typename CT>
		array_<N, CT> operator+(array_<N, CT> a, const array_<N, CT>& b)
		{
			return ((array_<N,CT>)a) += b;
		}

		template<size_t N, typename CT>
		array_<N, CT> operator-( array_<N, CT> a, const array_<N, CT>& b)
		{
			return ((array_<N, CT>)a) -= b;
		}

		template<size_t N, typename CT>
		CT dot(const array_<N, CT>& a, const array_<N, CT>& b)
		{
			return ((array_<N, CT>)a).dot(b);
		}

		template< typename CT>
		array_<3, CT>  cross(const array_<3, CT>& a, const array_<3, CT>& b)
		{
			array_<3, CT> tmp;
			return cross_product(a, b, tmp);
			
		}

		template< typename CT>
		CT  cross(const array_<2, CT>& a, const array_<2, CT>& b)
		{
			return a[0] * b[1] - a[1] * b[0];

		}

		template<size_t N, typename CT, typename VT>
		array_<N, CT> operator*(const array_<N, CT> a, VT val)
		{
			return ((array_<N, CT>)a) *= val;
		}

		template<size_t N, typename CT, typename VT>
		array_<N, CT> operator*(VT val,const array_<N, CT> a)
		{
			return ((array_<N, CT>)a) *= val;
		}

		
		template<size_t N, typename CT>
		array_<N, CT> operator-(const Point<N, CT>& a, const Point<N, CT>& b)
		{
			array_<N, CT> tmp(a);

			return ((array_<N, CT>)a) -= b;
		}







	}// end namespace GEO

}// end namespace MML