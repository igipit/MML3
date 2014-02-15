#ifndef _ORDERED_LIST_MAP_H_
#define _ORDERED_LIST_MAP_H_

#include"MML3-base_ordered_list.h"
//#define BOOST_DISABLE_THREADS
//#include "boost/pool/pool.hpp"
#include<map>
namespace MML3
{

	template<typename KT, typename VT>
	class ordered_list_map :public base_ordered_list<KT, VT>
	{
		typedef base_ordered_list<KT, VT>						mybase;
		
		
		
		
	
	public:
		typedef typename mybase::value_t						value_t;
		typedef typename mybase::key_t							key_t;
		typedef std::pair<const KT, VT>							key_value_pair_t;
		typedef std::allocator<key_value_pair_t>				allocator_t;
		typedef std::map<KT, VT, std::less<KT>, allocator_t>	container_t;
		typedef typename container_t::iterator					iterator;
		typedef typename container_t::const_iterator			const_iterator;

		
		//typedef boost::fast_pool_allocator<key_value_pair_t>	allocator_t;

		ordered_list_map(const allocator_t& alloc = allocator_t()) :cont_(alloc){};
		ordered_list_map(const ordered_list_map&) = default;
		ordered_list_map(ordered_list_map&& o) :cont_(std::move(o.cont_)){}
		~ordered_list_map() = default;

		//----------------------------------------------------------------
		//OVERRIDE OF PURE VIRTUAL FUNCTIONS
		virtual void destroy()
		{ 
			cont_.clear(); 
		}

		virtual value_t& 	value_at(key_t key)
		{
			return cont_[key];
		}

		virtual value_t*	find(key_t key)
		{
			auto it = cont_.find(key);
			if (it == cont_.end())
				return nullptr;
			else
				return &(it->second);
		}
		virtual const value_t*	find(key_t key)const
		{
			auto it = cont_.find(key);
			if (it == cont_.end())
				return nullptr;
			else
				return &(it->second);
		}

		/// returns the size of the list
		size_t			size()const
		{ 
			return cont_.size(); 
		}


		//-----------------------------------------------------
		// ADDED FUNCTIONS


		iterator				begin(){ return cont_.begin(); }
		const_iterator			begin()const{ return cont_.begin(); }
		const_iterator			end()const{ return cont_.end(); }
		
		static value_t&			value( iterator& it){ return it->second; }
		static const value_t&	value(const const_iterator& it){ return it->second; }
		static const key_t&		key( iterator& it){ return it->first; }
		static const key_t&		key(const const_iterator& it){ return it->first; }

	private:
		container_t cont_;
	};




}// end namespace
#endif