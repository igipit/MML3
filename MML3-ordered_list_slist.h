#ifndef _ORDERED_LIST_SLIST_H_
#define _ORDERED_LIST_SLIST_H_

#include"MML3-base_ordered_list.h"
//#define BOOST_DISABLE_THREADS
//#include "boost/pool/pool.hpp"
#include<forward_list>
namespace MML3
{

	template<typename KT, typename VT>
	class ordered_list_slist :public base_ordered_list<KT, VT>
	{
		typedef base_ordered_list<KT, VT>						mybase;
		
	
	public:
		typedef typename mybase::value_t				value_t;
		typedef typename mybase::key_t					key_t;
		typedef std::pair<const KT, VT>					key_value_pair_t;
		typedef std::allocator<key_value_pair_t>		allocator_t;
		typedef std::forward_list<key_value_pair_t, allocator_t>
														container_t;

		typedef typename container_t::iterator			iterator;
		typedef typename container_t::const_iterator	const_iterator;

		
		ordered_list_slist(const allocator_t& alloc = allocator_t()) :cont_(alloc){};
		ordered_list_slist(const ordered_list_slist&) = default;
		ordered_list_slist(ordered_list_slist&& o) :cont_(std::move(o.cont_)){}
		~ordered_list_slist() = default;

		//----------------------------------------------------------------
		//OVERRIDE OF PURE VIRTUAL FUNCTIONS
		virtual void destroy()
		{ 
			cont_.clear(); 
		}

		virtual value_t& 	value_at(key_t key)
		{
			auto end = cont_.end();
			auto it = cont_.begin();
			
			// se la lista e' vuota oppure se il primo elemento ha chiave
			// maggiore di key inserisco un nuovo elemento all'inizio
			if (it == end || it->first > key)
			{
				cont_.push_front(key_value_pair_t(key, value_t()));
				return cont_.front().second;
			}
			// adesso cerco l'elemento (*it) con chiave <= key
			for (auto it2 = std::next(it); it2 != end; ++it2)
			{
				if (it2->first > key)
					break;
				else
					++it;
			}
			
			if (it->first != key)
				it = cont_.insert_after(it, key_value_pair_t(key, value_t()));
			
			return it->second;
		}

		virtual value_t*	find(key_t key)
		{
			for (auto it = cont_.begin(), end = cont_.end(); it != end; ++it)
			if (it->first == key)
				return &(it->second);
			return nullptr;
		}
		virtual const value_t*	find(key_t key)const
		{
			for (auto it = cont_.begin(), end = cont_.end(); it != end; ++it)
			if (it->first == key)
				return &(it->second);
			
			return nullptr;
		}

		/// returns the size of the list
		size_t			size()const
		{ 
			size_t count = 0;
			for (auto it = cont_.begin(), end = cont_.end(); it != end; ++it)
				++count;
			return count;
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