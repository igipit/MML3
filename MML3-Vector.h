#pragma once
#include<memory>
#include<iostream>
#include<iomanip>
#include<cassert>
#include"MML3-config.h"

namespace MML3
{

	template<typename T>
	class Vector
	{
		enum: size_t{Base_=BASE::OFFSET};
	public:
		// CTORs
		Vector()= default;
		explicit Vector(size_t sz)							:p_(sz?new T[sz]:nullptr), sz_(sz){}
		Vector(const Vector& o)								:Vector(o.begin(),o.size()){}
		Vector(Vector&& o)									:p_(o.p_.release()), sz_(o.sz_){}
		Vector(const T* o, size_t sz)						:p_(sz?new T[sz]:nullptr), sz_(sz){ std::copy(o, o + sz, p_.get()); }
		Vector(const std::initializer_list<T>& s)			:Vector(s.begin(), s.size()){}
		~Vector()= default;
		Vector&		operator=(const Vector& o)				{ return assign(o.p_.get(), o.sz_);}
		Vector&		operator=(T val)						{ std::fill_n(p_.get(), sz_, val); return *this; }
		Vector&		operator=(const std::initializer_list<T>& s) { return assign(s.begin(), s.size()); }
		Vector&		assign(const T* o, size_t sz)			{ resize(sz); std::copy(o, o + sz, p_.get()); return *this; }
		void		swap(Vector& o)							{ std::swap(p_,o.p_); std::swap(sz_, o.sz_); }
		Vector&		resize(size_t sz)						{ if (sz_ != sz){Vector tmp(sz);swap(tmp);}  return *this; }
		void		free()									{ p_.reset(); sz_ = 0; }

		// ACCESSORS
		T*			begin()									{return p_.get(); }
		T*			end()									{return p_.get()+sz_; }
		const T*	begin()const							{return p_.get(); }
		const T*	end()const								{return p_.get() + sz_; }
		T&			operator[](size_t i)					{return p_[i];}
		const T&	operator[](size_t i)const				{return p_[i];}
		T&			operator()(size_t i);
		const T&	operator()(size_t i)const;


		Vector&		operator+=(const Vector& o);
		Vector&		operator-=(const Vector& o);
		Vector&		operator+();
		Vector		operator-()const;
		Vector&		operator*=(T val);
		Vector&		operator/=(T val);
		//returns a copy of the Matrix transformed by func
		Vector		apply(T func(T)) const;
		Vector&		transform(T func(T));


		// INFO
		size_t		size()const{ return sz_;}
		bool		print(std::ostream& os)const;

	private:
	    std::unique_ptr<T[]>  p_;

		size_t				sz_ = 0;

	};

	template<typename T>
	std::ostream& operator << (std::ostream& os, const Vector<T>& obj)
	{
		obj.print(os);
		return os;
	}

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
	inline T&	Vector<T>::operator()(size_t i)
	{
#ifdef MML3_TEST_INDEX_ON_ACCESS
		assert(i>=Base_ && i < sz_+Base_);
#endif
		return p_[i - Base_];
	}

	template<typename T>
	inline const T&	Vector<T>::operator()(size_t i)const
	{
#ifdef MML3_TEST_INDEX_ON_ACCESS
		assert(i >= Base_ && i < sz_ + Base_);
#endif
		return p_[i - Base_];
	}


	template<typename T>
	bool  Vector<T>::print(std::ostream& os)const
	{
		std::streamsize ssize = os.width();
		std::streamsize sprecision = os.precision();
		os << "[" << size() <<"]"<< std::endl;
		for (size_t r = 0; r != sz_; ++r)
			os << std::setw(ssize) << std::setprecision(sprecision) << p_[r] << ' ';
		os << std::endl;
		return os.good();
	}




	template<typename T>
	inline Vector<T>&	Vector<T>::operator+=(const Vector& o)
	{
		if (size() != o.size())
			throw std::runtime_error("Vector +=: size mismatch");
		T*			dest = p_.get();
		const T*	src = o.p_.get();
		for (size_t i = 0; i != sz_; ++i)
			dest[i] += src[i];
		return *this;
	}

	template<typename T>
	inline Vector<T>&	Vector<T>::operator-=(const Vector& o)
	{
		if (size() != o.size())
			throw std::runtime_error("Vector -=: size mismatch");
		T*			dest = p_.get();
		const T*	src = o.p_.get();
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
			tmp.p_[i] = -p_[i];
		return tmp;
	}

	template<typename T>
	inline Vector<T>&	Vector<T>::operator*=(T val)
	{
		for (size_t i = 0; i != sz_; ++i)
			p_[i]*=val;
		return *this;
	}

	template<typename T>
	inline Vector<T>&	Vector<T>::operator/=(T val)
	{
		for (size_t i = 0; i != sz_; ++i)
			p_[i] /= val;
		return *this;
	}

	template<typename T>
	inline Vector<T>	Vector<T>::apply(T func(T)) const
	{
		Vector tmp(sz_);
		for (size_t i = 0; i != sz_; ++i)
			tmp.p_[i] = func(p_[i]);
		return tmp;
	}

	template<typename T>
	inline Vector<T>&	Vector<T>::transform(T func(T))
	{

		for (size_t i = 0; i != sz_; ++i)
			p_[i] = func(p_[i]);
		return *this;
	}







}
