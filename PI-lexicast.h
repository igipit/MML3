
#pragma once
#include<type_traits>
#include<cstdint>
#include<string>
#include<stdexcept>
#include<sstream>
#include<typeinfo>
namespace PI
{

	//// converts anything to anything; throws if conversion is not possible
	//template<typename Target, typename Source>
	//Target lexicast(const Source& arg);

	

	// converts anything to a string
	template<typename T>
	std::string as_string(const T& val);


	// converts a string to a number of any type
	template<typename Number>
	Number as_number(const std::string&);
	// converts a string to a number
	template<typename T>
	T str2num(const std::string& arg){ return as_number<T>(arg); }


	////////////////////////////////////////////////////////////
	// IMPLEMENTATION
	template<typename A, typename B>
	struct arithmetic_A2B
	{
		static B cast(A a)
		{
			if (a > std::numeric_limits<B>::max() || a < std::numeric_limits<B>::min())
				throw std::out_of_range("PI::arithmetic_A2B::cast");
			return B(a);
		}
	};


	////////////////////////////////////////////////////////
	// STRING TO NUMBER CAST
		
	// throws std::out_of_range or std::invalid_argument
	template<typename Number>
	Number as_number(const std::string& arg)
	{

		if (std::is_same<Number,int>::value)
			return (Number)std::stoi(arg);

		if (std::is_same<Number, long>::value)
			return (Number)std::stol(arg);

		if (std::is_same<Number, unsigned long>::value)
			return (Number)std::stoul(arg);

		if (std::is_same<Number, long long>::value)
			return (Number)std::stoll(arg);

		if (std::is_same<Number, double>::value)
			return (Number)std::stod(arg);
		if (std::is_same<Number, float>::value)
			return (Number)std::stof(arg);
		
		if (! std::is_arithmetic<Number>::value)
			assert(false);

		if (std::is_integral<Number>::value)
		{
			long long value = std::stoll(arg);
			return (Number)arithmetic_A2B<long long, Number>::cast(value);
		}
		
		assert(false);
		return Number();
	}

	////  SPECIALIZED VERSIOS
	//// int
	//template<> 	
	//inline int				as_number<int>(const std::string& arg){	return std::stoi(arg);}
	//// long  
	//template<>
	//inline long				as_number<long>(const std::string& arg){ return std::stol(arg); }
	//// unsigned long
	//template<>
	//inline unsigned long	as_number<unsigned long>(const std::string& arg){ return std::stoul(arg); }
	//// double
	//template<>
	//inline double			as_number<double>(const std::string& arg){ return std::stod(arg); }
	//// float
	//template<> 	
	//inline float			as_number<float>(const std::string& arg){ return std::stof(arg); }



	

	// anything to string
	template<typename T>
	std::string as_string(const T& val)
	{
		if (std::is_arithmetic<T>::value)
			return std::to_string(val);
		std::ostringstream oss;
		oss << val;
		if (!oss)
			throw std::bad_cast("PI::as_string: cast failed");
		return oss.str();

	}


	// universal but generally slow lexicast (10-50 times slower than specialized versions)
	template<typename Target, typename Source>
	Target gen2gen_lexicast(const Source& arg)
	{
		std::istringstream iss(arg);
		Target val;
		iss >> val;
		if (iss.fail())
			throw std::range_error("");
		return val;
	}



	
	//// anything to anything
	//template<typename Target, typename Source>
	//Target lexicast(const Source& arg)
	//{
	//	return gen2gen_lexicast<Target>(arg);
	//}

	//// anything to string
	//template<typename Source>
	//std::string lexicast(const Source& arg){return as_string(arg);}


	//// string to any number
	//template<typename Target>
	//Target lexicast(const std::string& arg)
	//{
	//	static_assert(std::is_arithmetic<Target>::value,"");
	//	return as_number<Target>(arg);
	//}
	//// string to string
	//template<>
	//std::string lexicast<std::string>(const std::string& arg)
	//{
	//	return arg;
	//}


} // end namespace PI



