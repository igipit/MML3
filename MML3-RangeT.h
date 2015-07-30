
#pragma once

namespace MML3
{

	template<typename T=int>
	class Range
	{
		typedef T  							element_t;
	public:
		struct iterator{
			element_t v;
			element_t operator*()const{ return v; }
			operator element_t()const{ return v; }
			// pre increment op
			element_t operator++(){ return v += element_t(1); }
			bool       operator!=(const iterator& o)const{ return v != o.v; }
		};

		Range(element_t first, element_t last) :begin_(first), end_(last+element_t(1)){}
		Range() = delete;
		Range(const Range&) = delete;

		iterator begin()const{ return iterator{ begin_ }; }
		iterator end()const{ return iterator{ end_ }; }

		element_t		 max()const{ return end_ - 1;}
		element_t		 min()const{ return begin_ ;}
		size_t			 size()const{ return end_ - begin_;}
		element_t		operator[](size_t i)const{ return begin_ + i;}

	private:
		element_t begin_ = element_t(0), end_ = element_t(0);
	};
}// end namespace