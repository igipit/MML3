#ifndef _MML3_STATIC_ARRAY_H_
#define _MML3_STATIC_ARRAY_H_

#include"MML3-config.h"
#include "MML3-Math.h"
#include"MML3-Math.h"
#include<array>
#include<vector>
#include <cmath>

namespace MML3
{


template<typename T, size_t N1, size_t N2>
class StaticArray;


template<typename T, size_t N1>
using StaticVector=StaticArray<T,N1,1>;

template<typename T, size_t N1, size_t N2>
struct base_static_array_
{

	enum{ Base_ = 1 };

	typedef T	value_t;
	std::array<T, N1*N2>  p_;

	base_static_array_() = default;
	base_static_array_(const base_static_array_& o) = default;
	base_static_array_(base_static_array_&& o) :p_(std::move(o.p_)){}
	virtual ~base_static_array_() = default;

	T*			data(){ return p_.data(); }
	const T*	data()const{ return p_.data(); }
	size_t		size()const{ return N1*N2; }
	void		fill(const value_t& val){ p_.fill(val); }
};


template<typename T, size_t N1, size_t N2>
struct array_data_;


template <typename T, size_t M, size_t N>
struct array_data_ :public base_static_array_<T, M, N>
{
	typedef base_static_array_<T, M, N> base_;
	typedef T					value_t;

	array_data_() = default;
	array_data_(const array_data_& o) = default;
	array_data_(array_data_&& o) :base_(std::move(o)){}
	virtual ~array_data_() = default;


	T*			operator[](size_t i){ return &p_[i*N]; }
	const T*	operator[](size_t i)const{ return &p_[i*N]; }


	T&			operator()(size_t i, size_t j){ return at_0b(i - Base_, j - Base_); }
	const T&	operator()(size_t i, size_t j)const { return at_0b(i - Base_, j - Base_); }
	T&			at_0b(size_t i, size_t j){ return p_[i*N + j]; }
	const T&	at_0b(size_t i, size_t j)const { return p_[i*N + j]; }
};


template <typename T, size_t M>
struct array_data_<T, M, 1>:public base_static_array_<T, M, 1>
{
	typedef base_static_array_<T, M, 1> base_;
	typedef T					value_t;

	array_data_() = default;
	array_data_(const array_data_& o) = default;
	array_data_(array_data_&& o) :base_(std::move(o)){}
	virtual ~array_data_() = default;

	T&			operator[](size_t i){ return p_[i]; }
	const T&	operator[](size_t i)const{ return p_[i]; }


	T&			operator()(size_t i){ return p_[i - Base_]; }
	const T&	operator()(size_t i)const{ return p_[i - Base_]; }

	// second index ignored; provided for functions that access vectors as 
	// 1 column arrays
	T&			operator()(size_t i, size_t){ return p_[i - Base_]; }
	const T&	operator()(size_t i, size_t)const{ return p_[i - Base_]; }
	T&			at_0b(size_t i, size_t){ return p_[i]; }
	const T&	at_0b(size_t i, size_t)const { return p_[i]; }
};


template<typename T, size_t N1, size_t N2>
class StaticArray :public array_data_<T, N1, N2>
{
	typedef array_data_<T, N1, N2> base_class;
	typedef typename std::conditional<(N2>1), size_t, void > ::type idx2_t;
public:

	StaticArray() = default;
	StaticArray(const StaticArray&) = default;
	explicit StaticArray(T value){ fill(value); }
	StaticArray(const std::initializer_list<T>& s){ assign(s.begin(), s.size()); }
	StaticArray(const std::initializer_list<std::initializer_list<T>>& s){ assign(s); }
	StaticArray(const T* src, size_t src_sz){ assign(src, src_sz); }
	~StaticArray(){ base_class::~base_class(); }


	using base_class::data;
	//using base_class::end;
	using base_class::operator[];
	using base_class::operator();
	using base_class::at_0b;
	using base_class::value_t;

	StaticArray&	operator=(const StaticArray&) = default;
	StaticArray&	operator=(const std::initializer_list<T>& s){ return assign(s.begin(), s.size()); }
	StaticArray&	operator=(const std::initializer_list<std::initializer_list<T>>& s){ return assign(s); }
	StaticArray&	operator=(T value){ fill(value); return *this; }
	template<typename TT>
	StaticArray&	assign(const TT* src, size_t src_sz);
	StaticArray&	assign(const std::initializer_list<std::initializer_list<T>>& s);
	// provided for coherence with other array classes: throws if new sizes dont equal old sizes
	StaticArray&	resize(size_t N, size_t M){ if (N != N1 || M != N2) throw std::runtime_error("StaticArray not resizeable"); return *this; }

	size_t			size()const{ return N1*N2; }
	size_t          nrows()const{ return N1; }
	size_t          ncols()const{ return N2; }

	StaticArray&	operator+=(const StaticArray&);
	StaticArray&	operator-=(const StaticArray&);
	StaticArray&	operator*=(T);
	StaticArray&	operator/=(T);
	StaticArray&	operator+()									{ return *this; }
	StaticArray		operator-()const;

	bool			print(std::ostream& os)const;

	StaticArray&	set_identity();

};

// non member operators

template<typename T, size_t N1, size_t N2>
std::ostream& operator << (std::ostream& o, const StaticArray<T, N1, N2>& v)
{
	v.print(o);
	return o;
}

template<typename T, size_t N1, size_t N2>
inline StaticArray<T, N1, N2>   operator+(StaticArray<T, N1, N2> x, const StaticArray<T, N1, N2>& y)
{
	return x += y;
}

template<typename T, size_t N1, size_t N2>
inline StaticArray<T, N1, N2>   operator-(StaticArray<T, N1, N2> x, const StaticArray<T, N1, N2>& y)
{
	return x -= y;
}

template<typename T, size_t N1, size_t N2>
inline StaticArray<T, N1, N2>   operator*(StaticArray<T, N1, N2> x, T val)
{
	return x *= val;
}

template<typename T, size_t N1, size_t N2>
inline StaticArray<T, N1, N2>   operator*(T val, StaticArray<T, N1, N2> x)
{
	return x *= val;
}


template<typename T, size_t N1, size_t N2>
inline StaticArray<T, N1, N2>   operator/(const StaticArray<T, N1, N2> x, T val)
{
	return x /= val;
}



// array product

template<typename T, size_t M, size_t N, size_t K>
StaticArray<T, M, N> product(const StaticArray<T, M, K>& A, const StaticArray<T, K, N>& B)
{
	StaticArray<T, M, N> C;
	for (size_t m = 0; m != M; ++m)
	for (size_t n = 0; n != N; ++n)
	{
		T acc(0);
		for (size_t k = 0; k != K; ++k)
			acc += A.at_0b(m, k) * B.at_0b(k, n);
		C.at_0b(m, n) = acc;
	}
	return C;
}

// specialization for array - vector product
template<typename T, size_t M, size_t N>
StaticArray<T, M, 1> product(const StaticArray<T, M, N>& A, const StaticArray<T, N, 1>& B)
{
	StaticArray<T, M, 1> C;
	for (size_t m = 0; m != M; ++m)
	{
		T acc(0);
		for (size_t n = 0; n != N; ++n)
			acc += A[m][n] * B[n];
		C[m] = acc;
	}
	return C;
}


template<typename T, size_t M, size_t N>
StaticArray<T, N, M> transpose(const StaticArray<T, M, N>& A)
{
	StaticArray<T, N, M> tmp;
	for (size_t i = 0; i != N; ++i)
	for (size_t j = 0; j != M; ++j)
		tmp.at_0b(i, j) = A.at_0b(j, i);
	return tmp;
}


// array dot product
template<typename T, size_t M, size_t N>
inline T dot(const StaticArray<T, M, N>& A, const StaticArray<T, M, N >& B)
{

	T acc(0);
	const T* a = A.data();
	const T* b = B.data();
	for (size_t m = 0; m != M*N; ++m)
		acc += a[m] * b[m];
	return sqrt(acc);
}

// array dot product
template<typename T, size_t M, size_t N>
inline T norm2(const StaticArray<T, M, N>& A)
{
	return sqrt(dot(A, A));
}



// cross product in 3D space
template<typename T>
inline StaticArray<T, 3, 1> cross_product(const StaticArray<T, 3, 1>& A, const StaticArray<T, 3, 1 >& B)
{
	StaticArray<T, 3, 1> tmp;
	Math::vector_product_3D(A.data(), B.data(), tmp.data());
	return tmp;

}




// returns the determinant of a 2x2, 3x3, 4x4 array
template<typename T, size_t N>
inline T  det(const StaticArray<T, N, N>& A)
{
	return Math::Matrix<N, T>::det(A.data());
}

// invert a  2x2, 3x3, 4x4 array and returns the matrix determinanat

template<typename T, size_t N>
inline T  inv(StaticArray<T, N, N>& A)
{
	return Math::Matrix<3, T>::inv(A.data());
}








//////////////////////////////////////////////
// IMPLEMENTATION OF 2 DIMENSIONAL fs ARRAY
//////////////////////////////////////////////

template<typename T, size_t N1, size_t N2>
inline auto StaticArray<T, N1, N2>::set_identity()->StaticArray&
{
	assert(N1 == N2);
	*this = T(0);
	for (size_t i = 0; i != N1; ++i)
		at_0b(i, i) = T(1);
	return *this;
}


template<typename T, size_t N1, size_t N2>
inline auto StaticArray<T, N1, N2>::assign(const std::initializer_list<std::initializer_list<T>>& s)->StaticArray&
{
	size_t nr = std::min(N1, s.size());

	const std::initializer_list<T>* p_r = s.begin();
	for (size_t r = 0; r != nr; ++r)
	{

		size_t nc = std::min(N2, p_r->size());
		std::copy(p_r->begin(), p_r->begin() + nc, operator[](r));
		++p_r;
	}
	return *this;


}



template<typename T, size_t N1, size_t N2>
template<typename TT>
inline auto StaticArray<T, N1, N2>::assign(const TT* src, size_t src_sz)->StaticArray&
{
	// l'array attuale può essere più corto della sorgente e viceversa
	size_t n = std::min(size(), src_sz);
	T* p = data();
	for (size_t i = 0; i != n; ++i)
		p[i] = (T)src[i];
	return *this;
}


template<typename T, size_t N1, size_t N2>
bool StaticArray<T, N1, N2>::print(std::ostream& os)const
{
	std::streamsize ssize = os.width();
	std::streamsize sprecision = os.precision();
	os << "[" << nrows() << "," << ncols() << "] " << std::endl;;
	for (size_t r = 1; r <= N1; ++r)
	{
		for (size_t c = 1; c <= N2; ++c)
			os << std::setw(ssize) << std::setprecision(sprecision) << operator()(r, c) << ' ';
		os << std::endl;
	}
	return os.good();
}







template<typename T, size_t N1, size_t N2>
inline auto StaticArray<T, N1, N2>::operator+=(const StaticArray& o)->StaticArray&
{
	T* dest = data();
	const T* src = o.data();
	for (size_t i = 0; i != N1*N2; ++i)
		dest[i] += src[i];
	return *this;
}

template<typename T, size_t N1, size_t N2>
inline auto StaticArray<T, N1, N2>::operator-=(const StaticArray& o)->StaticArray&
{
	T* dest = data();
	const T* src = o.data();
	for (size_t i = 0; i != N1*N2; ++i)
		dest[i] -= src[i];
	return *this;
}

template<typename T, size_t N1, size_t N2>
inline auto StaticArray<T, N1, N2>::operator*=(T val)->StaticArray&
{
	T* dest = data();

	for (size_t i = 0; i != N1*N2; ++i)
		dest[i] *= val;
	return *this;
}

template<typename  T, size_t N1, size_t N2>
inline auto StaticArray<T, N1, N2>::operator/=(T val)->StaticArray&
{
	T* dest = data();

	for (size_t i = 0; i != N1*N2; ++i)
		dest[i] /= val;
	return *this;
}

template<typename T, size_t N1, size_t N2>
inline auto StaticArray<T, N1, N2>::operator-()const->StaticArray
{
	StaticArray tmp;
	T*			dest = tmp.data();
	const T*	src = data();
	for (size_t i = 0; i != N1*N2; ++i)
		dest[i] = -src[i];
	return tmp;

}

} // end namespace MML3
#endif
