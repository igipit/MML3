#pragma once
#include<type_traits>
#include<cstdint>



namespace MML3{


	template<typename T>
	struct type{
		

		static std::uint16_t id()
		{
			return std::uint16_t(
				sizeof(T)* 40 +
				int(std::is_arithmetic<T>::value) * 20 +
				int(std::is_floating_point<T>::value) * 10+
				int(std::is_integral<T>::value)*5+
				int(std::is_signed<T>::value));
		}
	};





} // end namespace