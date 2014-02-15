#pragma once
#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

#include"MML3-Vector.h"
#include<initializer_list>
#include<iomanip>
#include"MML3-config.h"

namespace MML3{

	class iSet
	{
	public:
		typedef index_type  				index_t;
		typedef Vector<index_t>				vector_t;
		typedef vector_t::iterator 			iterator;
		typedef vector_t::const_iterator	const_iterator;

		iSet() = default;
		iSet(const	iSet&  rhs) :data_(rhs.data_){};
		iSet(iSet&& rhs) :data_(std::move(rhs.data_)){};
		//iSet(index_type i) :data_(1, i ){}
		// creates a set of indexes {first, first+1,..., last}
		iSet(index_t first, index_t last) :data_(last - first + 1){ std::iota(data_.begin(), data_.end(), first); }
		template<typename T>
		// creates a sz-set of indexes initialized to  {p[0],p[1],...,p[sz-1]}
		iSet(const T* p, index_t sz) : data_(p, p + sz){}
		template<typename S>
		iSet(std::initializer_list<S> s) : data_(s.begin(), s.size()){}
		//allows the construction of a iSet from any iterable class (e.g  std::containers)
		template<typename A, typename B, template<typename, typename > class Container>
		iSet(const Container<A, B>& c) : data_(c.begin(), c.end()){}

		iSet& 				operator=(const iSet& ois) = default;
		iSet& 				operator=(iSet&& ois)			{ data_ = std::move(ois.data_);	return *this; }
		iSet& 				add(index_t i)					{ data_.push_back(i); return *this; }
		index_t&  			operator()(index_t i)			{ return data_[i-BASE::OFFSET]; }
		const index_t&   	operator()(index_t i)const		{ return data_[i-BASE::OFFSET]; }
		index_t&  			operator[](index_t i)			{ return data_[i]; }
		const index_t&   	operator[](index_t i)const		{ return data_[i]; }
		iSet&				fill(index_t v)					{ std::fill_n(data_.begin(), data_.size(), v); return *this;}




		size_t   			size()const						{ return data_.size(); }
		iSet&     			reserve(index_t cap)			{ data_.reserve(cap); return *this; }
		iSet&     			resize(index_t sz)				{ data_.resize(sz); return *this; }
		index_t				capacity()const					{ return data_.capacity(); }
		bool	  			print(std::ostream& os)const;
		iterator   			begin()							{ return data_.begin(); }
		const_iterator   	begin()const					{ return data_.begin(); }
		iterator   			end()							{ return data_.end(); }
		const_iterator   	end()const						{ return data_.end(); }

		index_t				max()const						{return *std::max_element(data_.begin(), data_.end());}
		index_t				min()const						{return *std::min_element(data_.begin(), data_.end()); }

		static iSet 		range(index_t first, index_t last) { return iSet(first, last); }
	private:
		vector_t data_;
	};

	inline std::ostream& operator<<(std::ostream& os, const iSet& V) { V.print(os);	return os; }




	class iRange
	{
		typedef index_type  							index_t;
	public:
		struct iterator{
			index_type v;
			index_type operator*()const{ return v; }
			operator index_type()const{ return v; }
			// pre increment op
			index_type operator++(){ return v+=index_type(1); }
			bool       operator!=(const iterator& o)const{ return v != o.v; }
		};

		iRange(index_t b, index_t endd) :begin_(b), end_(endd){}
		iterator begin()const{ return iterator{ begin_ }; }
		iterator end()const{ return iterator{ end_ }; }
		
	private:
		index_type begin_=0, end_=0;
	};



	

	/////////////////////////////////////////////////////////////////////////////////////////////////
	inline bool	iSet::print(std::ostream& os)const
	{
		std::streamsize ssize = os.width();
		os << size() << std::endl;
		for (index_t r = BASE::first(); r != BASE::end(size()); ++r)
			os << std::setw(ssize) << this->operator()(r) << ' ';
		return os.good();
	}


	
	


} // end namespace