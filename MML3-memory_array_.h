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

	public:
		// CTOR's  ----------------------
		memory_array_() = default;
		memory_array_(const memory_array_& rhs)				{ copy_from(rhs.p_, rhs.sz_); }
		memory_array_(memory_array_&& rhs)					{ swap(rhs); }
		explicit memory_array_(size_t sz) :p_(new T[sz]), sz_(sz){}
		template<typename S>
		memory_array_(const S* src, size_t sz)			{ copy_from(src, sz); }
		template<typename S>
		memory_array_(std::initializer_list<S> s)		{ copy_from(s.begin(), s.size()); }

		~memory_array_() { free(); }
		// METHODs ----------------------
		void		resize(size_t sz)			{ if (sz_ != sz){ memory_array_ tmp(sz); swap(tmp); } }
		void		free()						{ delete[]p_; p_ = nullptr; sz_ = 0; }
		memory_array_&		swap(memory_array_& rhs)			{ std::swap(p_, rhs.p_); std::swap(sz_, rhs.sz_); return *this; }
		size_t		size()const					{ return sz_; }

		T*			begin()						{return sz_?p_:nullptr;}
		const T*	begin()const				{return sz_?p_:nullptr;}
		T*			end()						{return sz_?p_+sz_:nullptr;}
		const T*	end()const					{return sz_?p_+sz_:nullptr;}

		T&			operator[](size_t i)		{ return p_[i]; }
		const T&	operator[](size_t i)const	{ return p_[i]; }
		memory_array_&     operator=(const memory_array_& rhs){ return copy_from(rhs.p_, rhs.sz_); }
		memory_array_&     operator=(memory_array_&& rhs)		{ free(); std::swap(p_, rhs.p_); std::swap(sz_, rhs.sz_); return *this; }
		template<typename S>
		memory_array_&     operator=(std::initializer_list<S> s)    { return copy_from(s.begin(), s.size()); }
		memory_array_&		fill(T val)					{ for (size_t i = 0; i != sz_; ++i) p_[i] = val; return *this; }

		memory_array_&     copy_from(const T* src, size_t sz)
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
