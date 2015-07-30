#pragma once
#include<memory>
#include<iostream>
#include<iomanip>
#include<cassert>
#include<vector>
#include"MML3-config.h"

namespace MML3
{

	/** A  Vector class for arithmetic types (like std::valarray). 
		It is derived from std::vector	and adds:
			1) component access by operator () with 1-based indexes
			2) the basic arithmetic operators +=,-=,/=,*=
			3) apply, transform and fill methos
		Component access with operator []  conforms with C++ rule and thus is 0-based
	*/
	template<typename T>
	class Vector:public std::vector<T>
	{
		enum: size_t{Base_=BASE::OFFSET};
		typedef std::vector<T>	base_;
	public:

		static_assert(std::is_arithmetic<T>::value, "MML3::Vector is only for arithmetic types");

		typedef typename base_::value_type		value_type;
		typedef typename base_::value_type		value_t;
		typedef typename base_::size_type		size_type;
		typedef typename base_::difference_type difference_type;
		typedef typename base_::pointer			pointer;
		typedef typename base_::const_pointer	const_pointer;
		typedef typename base_::reference		reference;
		typedef typename base_::const_reference const_reference;
		typedef typename base_::const_iterator	const_iterator;

		// CTORs
		Vector()= default;
		explicit Vector(size_t sz)								:base_(sz){}
		Vector(const Vector& o)									:base_(o){}
		Vector(Vector&& o)										:base_(o){}
		Vector(const T* o, size_t sz)							:base_(o, o + sz){}

		// ci pensa il costruttore di base_ a verificare se il tipo Iter e' un iteratore
		template<class Iter>
		Vector(Iter beg, Iter end) : base_(beg, end){}
	

		Vector(const std::initializer_list<T>& s)				:base_(s){}
		~Vector()= default;


		Vector&		resize(size_t sz){ base_::resize(sz); return *this;}


		// cast operator to std::vector<T>
		operator	base_(){ return *this; }
		Vector&		swap(Vector& o)								{ base_::swap(o); return *this; }
		Vector&		operator=(const Vector& o)					{ base_::operator=(o); return *this; }
		Vector&		operator=(T val)							{ return fill(val); }
		Vector&		operator=(const std::initializer_list<T>& s){ base_::operator=(s); return *this; }
		Vector&		assign(const T* o, size_t sz)				{ 
			resize(sz); 
			if(sz)
				std::copy(o, o + sz, data()); 
			return *this; 
		}
		
		
		Vector&		fill(T val)						{ std::fill_n(data(),size(), val); return *this; }

#ifdef MML3_TEST_INDEX_ON_ACCESS
		T&			operator[](size_t i)			{ return base_::at(i);}
		const T&	operator[](size_t i)const		{ return base_::at(i); }
		// ACCESSORS 1-BASED
		T&			operator()(size_t i)			{ return base_::at(i - Base_); }
		const T&	operator()(size_t i)const		{ return base_::at(i - Base_); }
#else
		T&			operator()(size_t i)			{ return base_::operator[](i - Base_); }
		const T&	operator()(size_t i)const		{ return base_::operator[](i - Base_); }
#endif

		// VALUE VECTOR OPERATORS
		Vector&		operator+=(const Vector& o);
		Vector&		operator-=(const Vector& o);
		Vector&		operator+();
		Vector		operator-()const;
		Vector&		operator*=(T val);
		Vector&		operator/=(T val);
		//returns a copy of the Matrix transformed by func
		Vector		apply(T func(T)) const;
		Vector&		transform(T func(T));
	};

	
	template<typename T>
	inline Vector<T> operator +(const Vector<T>& x, const Vector<T>& y) { return (Vector<T>(x) += y); }

	template<typename T>
	inline Vector<T> operator -(const Vector<T>& x, const Vector<T>& y) { return (Vector<T>(x) -= y); }

	template<typename T>
	inline Vector<T> operator *(const Vector<T>& x, T val) { return (Vector<T>(x) *= val); }

	template<typename T>
	inline Vector<T> operator *( T val, const Vector<T>& x) { return (Vector<T>(x) *= val); }

	template<typename T>
	inline Vector<T> operator /(const Vector<T>& x, T val) { return (Vector<T>(x) /= val); }





	//////////////////////////////////////////////////////////////////////////////////////////////
	///					IMPLEMENTATION
	//////////////////////////////////////////////////////////////////////////////////////////////
	


	template<typename T>
	inline Vector<T>&	Vector<T>::operator+=(const Vector& o)
	{
		if (size() != o.size())
			throw std::runtime_error("Vector +=: size mismatch");
		T*			dest = data();
		const T*	src = o.data();
		size_t sz_ = size();
		for (size_t i = 0; i != sz_; ++i)
			dest[i] += src[i];
		return *this;
	}

	template<typename T>
	inline Vector<T>&	Vector<T>::operator-=(const Vector& o)
	{
		if (size() != o.size())
			throw std::runtime_error("Vector -=: size mismatch");
		T*			dest = data();
		const T*	src = o.data();
		size_t sz_ = size();
		for (size_t i = 0; i != sz_; ++i)
			dest[i] -= src[i];
		return *this;
	}

	template<typename T>
	inline Vector<T>&	Vector<T>::operator+(){ return *this;}

	template<typename T>
	inline Vector<T>	Vector<T>::operator-()const
	{
		Vector tmp(sz_);
		for (size_t i = 0; i != sz_; ++i)
			tmp.at(i) = -at(i);
		return tmp;
	}

	template<typename T>
	inline Vector<T>&	Vector<T>::operator*=(T val)
	{
		size_t sz_ = size();
		for (size_t i = 0; i != sz_; ++i)
			at(i)*=val;
		return *this;
	}

	template<typename T>
	inline Vector<T>&	Vector<T>::operator/=(T val)
	{
		size_t sz_ = size();
		for (size_t i = 0; i != sz_; ++i)
			at(i) /= val;
		return *this;
	}

	template<typename T>
	inline Vector<T>	Vector<T>::apply(T func(T)) const
	{
		size_t sz_ = size();
		Vector tmp(sz_);
		for (size_t i = 0; i != sz_; ++i)
			tmp.at(i) = func(at(i));
		return tmp;
	}

	template<typename T>
	inline Vector<T>&	Vector<T>::transform(T func(T))
	{
		size_t sz_ = size();
		for (size_t i = 0; i != sz_; ++i)
			at(i) = func(at(i));
		return *this;
	}







}
