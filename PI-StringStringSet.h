#ifndef _PI_STRING_STRING_SET_H_
#define _PI_STRING_STRING_SET_H_


#include<string>
#include<map>
#include<algorithm>
#include<cctype>

/*inline std::string to_upper(const std::string& str)
{
std::string  s(str);
std::transform(s.begin(), s.end(), s.begin(), std::toupper);
}

*/



namespace PI
{



	// definisco un comparatote  di stringhe che ignora la differenza tra minuscolo e maiuscolo

	template <class CHT = char>
	struct str_iless;


	template <>
	struct str_iless<char>
	{
		typedef std::basic_string<char>	string_t;
		bool operator() (const string_t& x, const string_t& y) const
		{
			return (_stricmp(x.c_str(), y.c_str()) < 0);
		}
	};

	template <>
	struct str_iless<wchar_t>
	{
		typedef std::basic_string<wchar_t>	string_t;
		bool operator() (const string_t& x, const string_t& y) const
		{
			return (_wcsicmp(x.c_str(), y.c_str()) < 0);
		}
	};

	// comparazione tra stringhe ignorando il caso
	inline int str_icompare(const std::string& str1, const std::string& str2){ return _stricmp(str1.c_str(), str2.c_str()); }
	inline int str_icompare(const std::wstring& str1, const std::wstring& str2){ return _wcsicmp(str1.c_str(), str2.c_str()); }


	// mappa di stringhe stringhe che ignora il caso maiuscolo minuscolo
	template<typename CHT>
	using istring_istring_set = std::map<std::basic_string<CHT>, std::basic_string<CHT>, str_iless<CHT>>;

} // end namespace

#endif