#ifndef _ORDERED_LIST_MAP_H_
#define _ORDERED_LIST_MAP_H_

#include"MML3-base_ordered_list.h"
#include<map>
namespace MML3
{

	template<typename KT, typename VT>
	class ordered_list_map 
	{
		
	
	public:
		typedef VT												value_t;
		typedef KT												key_t;
		typedef std::pair<const KT, VT>							key_value_pair_t;
		typedef std::allocator<	key_value_pair_t>				allocator_t;
		typedef std::map<KT, VT, std::less<KT>, allocator_t>	container_t;
		typedef typename container_t::iterator					iterator;
		typedef typename container_t::const_iterator			const_iterator;

		
		//typedef boost::fast_pool_allocator<key_value_pair_t>	allocator_t;

		ordered_list_map(const allocator_t& alloc = allocator_t()) :cont_(alloc){};
		ordered_list_map(const ordered_list_map&) = default;
		ordered_list_map(ordered_list_map&& o) :cont_(std::move(o.cont_)){}
		~ordered_list_map() = default;

		//----------------------------------------------------------------
		// FUNDAMENTAL METHODS
		virtual void			destroy() { cont_.clear(); }
		virtual value_t& 		value_at(key_t key)	{return cont_[key];}
		virtual value_t*		find(key_t key)

		{
			auto it = cont_.find(key);
			return (it == cont_.end())?nullptr : &(it->second);
		}

		virtual const value_t*	find(key_t key)const
		{
			auto it = cont_.find(key);
			return (it == cont_.end()) ? nullptr : &(it->second);
		}

		/// returns the size of the list
		size_t					size()const	{ return cont_.size(); }


		iterator				begin(){ return cont_.begin(); }
		const_iterator			begin()const{ return cont_.begin(); }
		iterator				end(){ return cont_.end(); }
		const_iterator			end()const{ return cont_.end(); }


		//---------------------------------------------------------------------------------------------
		// implementation based on fundaental methods
		value_t &          put(const key_t& key, const value_t& value)	{
			return value_at(key) = value;
		}

		value_t &          sum(const key_t& key, const value_t& value)	{
			return value_at(key) += value;
		}

		// N.B. since template method are not allowed to be virtual here are some specializations
		template<typename VT1, typename IT1>
		void				put(const IT1* key, const VT1* value, size_t sz)	{
			for (size_t i = 0; i != sz; ++i)  value_at(key_t(key[i])) = value_t(value[i]);
		}

		template<typename VT1, typename IT1>
		void				sum(const IT1* key, const VT1* value, size_t sz)	{
			for (size_t i = 0; i != sz; ++i)  value_at(key_t(key[i])) += value_t(value[i]);
		}



		//-----------------------------------------------------
		// ADDED FUNCTIONALITY

		
		// converters from iterator to key and value
		static value_t&			value( iterator& it){ return it->second; }
		static const value_t&	value(const const_iterator& it){ return it->second; }
		static const key_t&		key( iterator& it){ return it->first; }
		static const key_t&		key(const const_iterator& it){ return it->first; }

	private:
		container_t cont_;
	};




}// end namespace
#endif