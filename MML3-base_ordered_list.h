#ifndef _BASE_ORDERED_LIST_H_
#define _BASE_ORDERED_LIST_H_
#include<cstdint>
namespace MML3
{
	template<typename KT, typename VT>
	class base_ordered_list
	{
	public:
		typedef KT			key_t;
		typedef VT			value_t;
		
		base_ordered_list() = default;
		base_ordered_list(const base_ordered_list& ) = default;
		virtual ~base_ordered_list(){};

		///-------------------------------------------------------
		/// pure virtual functions
		///-------------------------------------------------------

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
		// distrugge il contenuto
		virtual void			destroy() = 0;

		
	public:
		//---------------------------------------------------------------------------------------------
		// implementation based on virtual functions (YOU CAN OVERRIDE but it is not usually necessary)
		virtual value_t &          put(const key_t& key, const value_t& value)	{
			return value_at(key) = value;										}

		virtual value_t &          sum(const key_t& key, const value_t& value)	{
			return value_at(key)+=value;										}

		//---------------------------------------------------------------------------------------------
		// implementation based on virtual functions (YOU SHOULD  OVERRIDE to gain efficiency)

		// N.B. since template method are not allowed to be virtual here are some specializations
		virtual void				put(const std::int32_t* key, const double* value, size_t sz)	{
					for (size_t i = 0; i != sz; ++i)  value_at(key_t(key[i])) = value_t(value[i]);}

		virtual void				put(const std::int64_t* key, const double* value, size_t sz)	{
			for (size_t i = 0; i != sz; ++i)	value_at(key_t(key[i])) = value_t(value[i]);		}

		virtual void				put(const std::uint32_t* key, const double* value, size_t sz)	{
			for (size_t i = 0; i != sz; ++i)		value_at(key_t(key[i])) = value_t(value[i]);	}

		virtual void				put(const std::uint64_t* key, const double* value, size_t sz)	{
			for (size_t i = 0; i != sz; ++i)		value_at(key_t(key[i])) = value_t(value[i]);	}


		virtual void				sum(const std::int32_t * key, const double* value, size_t sz)	{
			for (size_t i = 0; i != sz; ++i)		value_at(key_t(key[i])) += value_t(value[i]);	}

		virtual void				sum(const std::uint32_t * key, const double* value, size_t sz)	{
			for (size_t i = 0; i != sz; ++i)		value_at(key_t(key[i])) += value_t(value[i]);	}

		virtual void				sum(const std::int64_t * key, const double* value, size_t sz)	{
			for (size_t i = 0; i != sz; ++i)		value_at(key_t(key[i])) += value_t(value[i]);	}

		virtual void				sum(const std::uint64_t * key, const double* value, size_t sz)	{
			for (size_t i = 0; i != sz; ++i)		value_at(key_t(key[i])) += value_t(value[i]);	}

	};





}// end namespace
#endif