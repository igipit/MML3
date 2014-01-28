#ifndef _MML3_STATIC_ARRAY_H_
#define _MML3_STATIC_ARRAY_H_

#include"MML3-config.h"
#include<array>
//#include<vector>
#include <cmath>

namespace MML3
{

template<class T, size_t N1, size_t...  Sizes>
class st_array;



template <typename T, size_t M, size_t N>
struct array_data_
{
	enum{ Base_ = BASE::OFFSET };
	T p_[M][N];
	T*			begin(){ return p_[0];}
	const T*	begin()const{ return p_[0]; }

	T*			end(){ return p_[M] + N; }
	const T*	end()const{ return p_[M] + N; }

	T*			operator[](size_t i){ return p_[i]; }
	const T*	operator[](size_t i)const{ return p_[i]; }
	T&			operator()(size_t i, size_t j){ return p_[i - Base_][j - Base_]; }
	const T&	operator()(size_t i, size_t j)const { return p_[i - Base_][j - Base_]; }

	T&			at_0b(size_t i, size_t j){ return p_[i][j]; }
	const T&	at_0b(size_t i, size_t j)const { return p_[i][j];}
};


template <typename T, size_t M>
struct array_data_<T,M,1>
{
	enum{ Base_ = BASE::OFFSET };
	T p_[M];
	T*			begin(){ return p_; }
	const T*	begin()const{ return p_; }
	T*			end(){ return p_ + M; }
	const T*	end()const{ return p_ + M; }
	T&			operator[](size_t i){ return p_[i]; }
	const T&	operator[](size_t i)const{ return p_[i]; }
	T&			operator()(size_t i){ return p_[i - Base_]; }
	const T&	operator()(size_t i)const{ return p_[i-Base_]; }
	// second index ignored
	T&			operator()(size_t i, size_t j){ return p_[i - Base_]; }
	const T&	operator()(size_t i, size_t j)const{ return p_[i - Base_]; }

	T&			at_0b(size_t i, size_t){ return p_[i];}
	const T&	at_0b(size_t i, size_t)const { return p_[i];}
};


//---------------------------------------------------
// static array class
//---------------------------------------------------
template<class T, size_t N1, size_t N2>
class st_array:public array_data_<T,N1,N2>
{
	typedef array_data_<T, N1, N2> base_class;
	typedef typename std::conditional<(N2>1), size_t, void > ::type idx2_t;
	
public:

	

	st_array() = default;
	st_array(const st_array&) = default;
	st_array(T value){ std::fill_n(begin(), size(), value); }
	st_array(const std::initializer_list<T>& s){ assign(s.begin(),s.size()); }
	st_array(const std::initializer_list<std::initializer_list<T>>& s){ assign(s); }
	st_array(const T* src, size_t src_sz){ assign(src, src_sz); }
	~st_array() = default;



	using base_class::begin;

	using base_class::end;

	using base_class::operator[];

	using base_class::operator();

	
	using base_class::at_0b;

	

	st_array&	operator=(const st_array&) = default;
	st_array&	operator=(const std::initializer_list<T>& s){ return assign(s.begin(), s.size()); }
	st_array&	operator=(const std::initializer_list<std::initializer_list<T>>& s){ return assign(s); }
	st_array&	operator=(T value){ return fill(value);}
	
	st_array&	assign(const T* src, size_t src_sz);
	st_array&	assign(const std::initializer_list<std::initializer_list<T>>& s);
	st_array&	fill(T value){ std::fill_n(begin(), size(), value); return *this; }

	size_t		size()const{ return N1*N2;}
	size_t      nrows()const{ return N1; }
	size_t      ncols()const{ return N2; }


	st_array&	operator+=(const st_array&);
	st_array&	operator-=(const st_array&);
	st_array&	operator*=(T);
	st_array&	operator/=(T);
	st_array&	operator+()									{ return *this; }
	st_array	operator-()const;
	
	bool		print(std::ostream& os)const;

	st_array&	set_identity()
	{ 
		assert(N1 == N2);  
		*this = T(0);
		for (size_t i = 0; i != N1; ++i)
			at_0b(i, i) = T(1);
		return *this;
	}

};

// non member operators

template<class T, size_t N1, size_t N2>
std::ostream& operator << (std::ostream& o, const st_array<T, N1, N2>& v)
{
	v.print(o);
	return o;
}

template<class T, size_t N1, size_t N2>
inline st_array<T, N1, N2>   operator+(st_array<T, N1, N2> x, const st_array<T, N1, N2>& y)
{
	return x += y;
}

template<class T, size_t N1, size_t N2>
inline st_array<T, N1, N2>   operator-(st_array<T, N1, N2> x, const st_array<T, N1, N2>& y)
{
	return x -= y;
}

template<class T, size_t N1, size_t N2>
inline st_array<T, N1, N2>   operator*(st_array<T, N1, N2> x, T val)
{
	return x *= val;
}

template<class T, size_t N1, size_t N2>
inline st_array<T, N1, N2>   operator*(T val, st_array<T, N1, N2> x)
{
	return x *= val;
}


template<class T, size_t N1, size_t N2>
inline st_array<T, N1, N2>   operator/(const st_array<T, N1, N2> x, T val)
{
	return x /= val;
}



// array product

template<typename T, size_t M1, size_t M2, size_t N1, size_t N2>
st_array<T, 1, 1> product(const st_array<T, M1, M2>, const st_array<T, N1, N2>& B)
{
	static_assert(false, "Not allowed: array size mismatch ");
}


template<typename T, size_t M, size_t N, size_t K>
st_array<T, M, N> product(const st_array<T, M, K>& A, const st_array<T, K, N>& B)
{
	st_array<T, M, N> C;
	for (size_t m = 0; m != M; ++m)
	for (size_t n = 0; n != N; ++n)
	{
		T acc(0);
		for (size_t k = 0; k != K; ++k)
			acc += A.at_0b(m,k) * B.at_0b(k,n);
		C.at_0b(m,n) = acc;
	}
	return C;
}

// specialization for array - vector product
template<typename T, size_t M, size_t N>
st_array<T,M,1> product(const st_array<T, M, N>& A, const st_array<T,N,1>& B)
{
	st_array<T, M,1> C;
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
st_array<T, N, M> transpose(const st_array<T, M, N>& A)
{
	st_array<T, N, M> tmp;
	for (size_t i = 0; i != N; ++i)
	for (size_t j = 0; j != M; ++j)
		tmp.at_0b(j,i) = A.at_0b(i,j);
	return tmp;
}


// array product
template<typename T, size_t M, size_t N>
T dot(const st_array<T, M, N>& A, const st_array<T, M, N >& B)
{

	T acc(0);
	const T* v = begin();
	for (size_t m = 0; m != M*N; ++m)
		acc += v[m] * v[m];
	return sqrt(acc);
}





//////////////////////////////////////////////
// IMPLEMENTATION OF 2 DIMENSIONAL fs ARRAY
//////////////////////////////////////////////



template<class T, size_t N1, size_t N2>
inline auto st_array<T, N1, N2>::assign(const std::initializer_list<std::initializer_list<T>>& s)->st_array&
{ 
	size_t nr = std::min(N1, s.size());

	const std::initializer_list<T>* p_r = s.begin();
	for (size_t r = 0; r != nr; ++r)
	{

		size_t nc = std::min(N2, p_r->size());
		std::copy(p_r->begin(), p_r->begin() + nc, p_[r]);
		++p_r;
	}
		return *this;
	
	
}



template<class T, size_t N1, size_t N2>
inline auto st_array<T, N1,N2>::assign(const T* src, size_t src_sz)->st_array&
{
	std::copy(src, src + std::max(size(), src_sz), begin());
	return *this;
}


template<class T, size_t N1, size_t N2>
bool st_array<T, N1,N2>::print(std::ostream& os)const
{
	std::streamsize ssize = os.width();
	std::streamsize sprecision = os.precision();
	os << "[" << nrows() << "," << ncols() << "] " << std::endl;;
	for (size_t r = 1; r <= N1; ++r)
	{
		for (size_t c = 1; c <= N2; ++c)
			os << std::setw(ssize) << std::setprecision(sprecision) << operator()(r,c) << ' ';
		os << std::endl;
	}
	return os.good();
}







template<class T, size_t N1, size_t N2>
inline auto st_array<T, N1, N2>::operator+=(const st_array& o)->st_array&
{
	T* dest			= begin();
	const T* src	= o.begin();
	for (size_t i = 0; i != N1*N2; ++i)
		dest[i] += src[i];
	return *this;
}

template<class T, size_t N1, size_t N2>
inline auto st_array<T, N1, N2>::operator-=(const st_array& o)->st_array&
{
	T* dest = begin();
	const T* src = o.begin();
	for (size_t i = 0; i != N1*N2; ++i)
		dest[i] -= src[i];
	return *this;
}

template<typename T, size_t N1, size_t N2>
inline auto st_array<T, N1, N2>::operator*=(T val)->st_array&
{
	T* dest = begin();
	
	for (size_t i = 0; i != N1*N2; ++i)
		dest[i] *= val;
	return *this;
}

template<typename  T, size_t N1, size_t N2>
inline auto st_array<T, N1, N2>::operator/=(T val)->st_array&
{
	T* dest = begin();

	for (size_t i = 0; i != N1*N2; ++i)
		dest[i] /= val;
	return *this;
}

template<class T, size_t N1, size_t N2>
inline auto st_array<T, N1, N2>::operator-()const->st_array
{
	st_array tmp;
	T*			dest	= tmp.begin();
	const T*	src		= begin();
	for (size_t i = 0; i != N1*N2; ++i)
		dest[i] = -src[i];
	return tmp;

}

} // end namespace MML3
#endif