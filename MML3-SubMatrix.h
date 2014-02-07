#pragma once

#include"MML3-config.h"
#include"MML3-iSet.h"
#include<stdexcept>

namespace MML3
{

template<typename T, typename MP, typename MS, typename MO>
class Matrix;
template<typename T, typename MAT>
class const_SubMatrix;

template<typename T, typename MAT>
class SubMatrix
{
	typedef MAT       Matrix_t;
	template<typename TT, typename MP, typename MS, typename MO>
	friend class Matrix;

	friend class const_SubMatrix<T,MAT>;
public:
	typedef index_type          index_t;
	typedef T                   value_t;

	SubMatrix(const SubMatrix&) = default;
	SubMatrix(SubMatrix&& rhs) :A_(rhs.A_), r_(std::move(rhs.r_)), c_(std::move(rhs.c_)){};
	// access operators for vectors
	T&			operator()(index_t r, index_t c)		{ return A_(r_(r), c_(c)); }
	const T&	operator()(index_t r, index_t c)const	{ return A_(r_(r), c_(c)); }
	T&			operator()(index_t r)					{ return A_(r_(r)); }
	const T&	operator()(index_t r)const				{ return A_(r_(r)); }
	index_t		nrows()const							{ return r_.size(); }
	index_t		ncols()const							{ return c_.size(); }
	index_t		size()const								{ return r_.size()*c_.size(); }


	SubMatrix&	fill( T val);
	SubMatrix&	operator=(const SubMatrix& rhs) = default;
	SubMatrix&	operator=(T val) {return fill(val);}
	
	template<typename MAT2>
	SubMatrix&	operator=(const MAT2& rhs);

	template<typename MAT2>
	SubMatrix&	operator+=(const MAT2& rhs);

	template<typename MAT2>
	SubMatrix&	operator-=(const MAT2& rhs);

	
	SubMatrix&	operator*=(const value_t& v);



	bool  print(std::ostream& os)const;

protected:
	SubMatrix(Matrix_t& A, const iSet& r, const iSet& c) :A_(A), r_(r), c_(c){ }
	SubMatrix(Matrix_t& A, iSet&& r, const iSet& c) :A_(A), r_(std::move(r)), c_(c){ }
	SubMatrix(Matrix_t& A, const iSet& r, iSet&& c) :A_(A), r_(r), c_(std::move(c)){ }
	SubMatrix(Matrix_t& A, iSet&& r, iSet&& c) :A_(A), r_(std::move(r)), c_(std::move(c)){ }
	Matrix_t&	A_;
	iSet	r_, c_;
};

template<typename T, typename MAT>
std::ostream& operator << (std::ostream& os, const SubMatrix<T, MAT>& obj)
{
	obj.print(os);
	return os;
}





template<typename T, class MAT>
class const_SubMatrix
{
	template<typename TT, typename MP, typename MS, typename MO>
	friend class Matrix;

	friend class SubMatrix<T, MAT>;
public:
	typedef index_type          index_t;
	typedef T                   value_t;
	typedef MAT					Matrix_t;
	typedef SubMatrix<T, MAT>   SubMatrix_t;
	//enum  SHAPE :bool { LT = MP::SHAPE_LT, GE = MP::SHAPE_GE };

	const_SubMatrix(const const_SubMatrix&) = default;
	const_SubMatrix(const SubMatrix_t& rhs) : A_(rhs.A_), r_(rhs.r_), c_(rhs.c_){};
	const_SubMatrix(SubMatrix_t&& rhs) : A_(rhs.A_), r_(std::move(rhs.r_)), c_(std::move(rhs.c_)){};
	// access operators for vectors
	const T&	operator()(index_t r, index_t c)const	{ return A_(r_(r), c_(c)); }
	const T&	operator()(index_t r)const				{ return A_(r_(r)); }
	index_t		nrows()const							{ return r_.size(); }
	index_t		ncols()const							{ return c_.size(); }
	index_t		size()const								{ return r_.size()*c_.size(); }


	bool  print(std::ostream& os)const;

protected:
	const_SubMatrix(const Matrix_t& A, const iSet& r, const iSet& c) :A_(A), r_(r), c_(c){ }
	const_SubMatrix(const Matrix_t& A, iSet&& r, const iSet& c) :A_(A), r_(std::move(r)), c_(c){ }
	const_SubMatrix(const Matrix_t& A, const iSet& r, iSet&& c) :A_(A), r_(r), c_(std::move(c)){ }
	const_SubMatrix(const Matrix_t& A, iSet&& r, iSet&& c) :A_(A), r_(std::move(r)), c_(std::move(c)){ }
private:
	const Matrix_t&	A_;
	iSet	r_, c_;
};

template<typename T, typename MAT>
std::ostream& operator << (std::ostream& os, const const_SubMatrix<T, MAT>& obj)
{
	obj.print(os);
	return os;
}





//print a matrix in text format (human readable),
template<typename T, typename MAT>
bool  SubMatrix<T,MAT>::print(std::ostream& os)const
{
	std::streamsize ssize = os.width();
	std::streamsize sprecision = os.precision();
	os << nrows() << " * " << ncols() << std::endl;

	for (auto r:r_)
	{
		for (auto c:c_)
		{
				os << std::setw(ssize) << std::setprecision(sprecision) << A_(r, c) << ' ';
		}
		os << std::endl;
	}
	return os.good();
}


//print a matrix in text format (human readable),
template<typename T, typename MAT>
bool  const_SubMatrix<T, MAT>::print(std::ostream& os)const
{
	std::streamsize ssize = os.width();
	std::streamsize sprecision = os.precision();
	os << nrows() << " * " << ncols() << std::endl;

	for (auto r : r_)
	{
		for (auto c : c_)
		{
			os << std::setw(ssize) << std::setprecision(sprecision) << A_(r, c) << ' ';
		}
		os << std::endl;
	}
	return os.good();
}

template<typename T, typename MAT>
auto SubMatrix<T, MAT>::fill(const T val)->SubMatrix&
{
	for (auto r : r_)
	for (auto c : c_)
		A_(r, c) = val;
	return *this;
}


template<typename T, typename MAT>
template<typename MAT2>
auto SubMatrix<T, MAT>::operator=(const MAT2& rhs)->SubMatrix&
{
	if (nrows() != rhs.nrows() || ncols() != rhs.ncols())
		throw std::range_error("SubMatrix<T, MAT>::operator= : range error");

	size_t nr = nrows(),nc=ncols();
	for (size_t i = 0; i != nr; ++i)
	for (size_t j = 0; j != nc; ++j)
		A_(r_[i], c_[j]) = rhs(i+1, j+1);
	return *this;
}


template<typename T, typename MAT>
template<typename MAT2>
auto SubMatrix<T, MAT>::operator+=(const MAT2& rhs)->SubMatrix&
{
	if (nrows() != rhs.nrows() || ncols() != rhs.ncols())
		throw std::range_error("SubMatrix<T, MAT>::operator+=: range error");

	size_t nr = nrows(), nc = ncols();
	for (size_t i = 0; i != nr; ++i)
	for (size_t j = 0; j != nc; ++j)
		A_(r_[i], c_[j])+= rhs(i + 1, j + 1);
	return *this;
}

template<typename T, typename MAT>
template<typename MAT2>
auto SubMatrix<T, MAT>::operator-=(const MAT2& rhs)->SubMatrix&
{
	if (nrows() != rhs.nrows() || ncols() != rhs.ncols())
		throw std::range_error("SubMatrix<T, MAT>::operator-=: range error");

	size_t nr = nrows(), nc = ncols();
	for (size_t i = 0; i != nr; ++i)
	for (size_t j = 0; j != nc; ++j)
		A_(r_[i], c_[j]) -= rhs(i + 1, j + 1);
	return *this;
}

template<typename T, typename MAT>
auto SubMatrix<T, MAT>::operator*=(const value_t& value)->SubMatrix&
{
	size_t nr = nrows(), nc = ncols();
	for (size_t i = 0; i != nr; ++i)
	for (size_t j = 0; j != nc; ++j)
		A_(r_[i], c_[j]) *= value;
	return *this;
}





} // end namespace
