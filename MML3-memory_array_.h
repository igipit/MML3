#pragma once
#include<memory>
#include<cstring>
//#include<type_traits>
#include<initializer_list>


namespace MML3{

	template<typename T>
	class memory_array_
	{
		//static_assert(std::is_trivially_copy_constructible<T>::value, "array_ is suited only for arithmetic types");
		typedef memory_array_ this_t;
	public:
		// CTOR's  ----------------------
		this_t() = default;
		this_t(const this_t& rhs)				{ copy_from(rhs.p_, rhs.sz_); }
		this_t(this_t&& rhs)					{ swap(rhs); }
		explicit this_t(size_t sz) :p_(new T[sz]), sz_(sz){}
		template<typename S>
		this_t(const S* src, size_t sz)			{ copy_from(src, sz); }
		template<typename S>
		this_t(std::initializer_list<S> s)		{ copy_from(s.begin(), s.size()); }

		~memory_array_() { free(); }
		// METHODs ----------------------
		void		resize(size_t sz)			{ if (sz_ != sz){ this_t tmp(sz); swap(tmp); } }
		void		free()						{ delete[]p_; p_ = nullptr; sz_ = 0; }
		this_t&		swap(this_t& rhs)			{ std::swap(p_, rhs.p_); std::swap(sz_, rhs.sz_); return *this; }
		size_t		size()const					{ return sz_; }

		T*			begin()						{return sz_?p_:nullptr;}
		const T*	begin()const				{return sz_?p_:nullptr;}
		T*			end()						{return sz_?p_+sz_:nullptr;}
		const T*	end()const					{return sz_?p_+sz_:nullptr;}

		T&			operator[](size_t i)		{ return p_[i]; }
		const T&	operator[](size_t i)const	{ return p_[i]; }
		this_t&     operator=(const this_t& rhs){ return copy_from(rhs.p_, rhs.sz_); }
		this_t&     operator=(this_t&& rhs)		{ free(); std::swap(p_, rhs.p_); std::swap(sz_, rhs.sz_); return *this; }
		template<typename S>
		this_t&     operator=(std::initializer_list<S> s)    { return copy_from(s.begin(), s.size()); }
		this_t&		fill(T val)					{ for (size_t i = 0; i != sz_; ++i) p_[i] = val; return *this; }

		this_t&     copy_from(const T* src, size_t sz)
		{
			resize(sz);
			std::memcpy( (void*)p_, (const void*)src, sizeof(T)*sz);
			return *this;
		}

	private:
		T		*p_  =nullptr;
		size_t	sz_ = 0;
	};


}// end namespace
