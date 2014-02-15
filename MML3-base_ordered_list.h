#ifndef _BASE_ORDERED_LIST_H_
#define _BASE_ORDERED_LIST_H_
namespace MML3
{
	template<typename KT, typename VT>
	class base_ordered_list
	{
	public:
		typedef KT	key_t;
		typedef VT	value_t;

		base_ordered_list() = default;
		base_ordered_list(const base_ordered_list& ) = default;
		virtual ~base_ordered_list(){};

		///-------------------------------------------------------
		/// pure virtual functions
		/// returns a reference to the value associated to key
		/// if the element does not exists it is created and intitialized with the default constructor
		virtual value_t& 		value_at(key_t key)=0;
		/// search for the element associated to the prescribed key
		///@return:  nullptr if the element does not exists
		///			 pointer to value else 	
		virtual value_t*		find(key_t key)=0;
		virtual const value_t*	find(key_t key)const = 0;
		/// returns the size of the list
		virtual size_t			size()const=0;

		virtual void destroy() = 0;

	public:
		// implementation based on virtual functions (you can override for efficiency)
		virtual value_t &          put(const key_t& key, const value_t& value)
		{
			return value_at(key) = value;
		}

		virtual value_t &          sum(const key_t& key, const value_t& value)
		{
			return value_at(key)+=value;
		}

		virtual void				put(const KT* key, const VT* value, size_t sz)
		{
			for (size_t i = 0; i != sz; ++i)
				value_at(key[i]) = value[i];
		}

		virtual void				sum(const KT* key, const VT* value, size_t sz)
		{
			for (size_t i = 0; i != sz; ++i)
				value_at(key[i]) += value[i];
		}

	};





}// end namespace
#endif