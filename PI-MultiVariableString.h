#ifndef _PI_MULTI_VARIABLE_STRING_H_
#define _PI_MULTI_VARIABLE_STRING_H_

#include<string>
#include<map>
#include <regex>
#include<tchar.h>
#include<algorithm>

#include"PI-StringStringSet.h"

namespace PI
{

	// definisco un comparatore di stringhe che ignora il caso minuscolo o maiuscolo





	

		


		///classe multi variable string
		///gestisce il formato MVS




	template<typename CHT>
	class basic_mvs;

	template<>
	class basic_mvs<char>
	{
	public:
		typedef char						char_t;
		typedef std::basic_string<char_t>	string_t;
		
		template<typename T>
		static string_t to_string(const T& v){ return std::to_string(v); }
		
		
		static const string_t& isep()
		{
			static string_t c("#");
			return c;
		}

		static const string_t& equal()
		{
			static string_t c("=");
			return c;
		}
		static const string_t& comma()
		{
			static string_t c(",");
			return c;
		}


		static const string_t& open_par()
		{
			static string_t c("(");
			return c;
		}
		static const string_t& close_par()
		{
			static string_t c(")");
			return c;
		}


	};

	template<>
	class basic_mvs<wchar_t>
	{
	public:
		typedef wchar_t					char_t;
		typedef std::basic_string<char_t>		string_t;

		template<typename T>
		static string_t to_string(const T& v){ return std::to_wstring(v); }

		static const string_t& isep()
		{
			static string_t c(L"#");
			return c;
		}

		static const string_t& equal()
		{
			static string_t c(L"=");
			return c;
		}
		static const string_t& comma()
		{
			static string_t c(L", ");
			return c;
		}

		static const string_t& open_par()
		{
			static string_t c(L"(");
			return c;
		}
		static const string_t& close_par()
		{
			static string_t c(L")");
			return c;
		}


	};


// NB: le stringhe vengono tutte convertite in small case
template<typename CHT>
class multi_variable_string :public basic_mvs<CHT>
{
	typedef basic_mvs<CHT>		base;
public:
	typedef typename base::char_t		char_t;
	typedef typename base::string_t     string_t;
	typedef istring_istring_set<char_t>	ssmap;
	typedef ssmap						var_set;


	/// genera una variabile indicizzata nel formato "var(i)"
	///@ param var			il nome della variabile
	///@ param i			l'indice 
	static string_t make_indexed_var(const char_t* str, int i)
	{
		string_t rv;
		rv.reserve(strlen(str) + 6);
		return rv.append(str).append(open_par()).append(to_string(i)).append(close_par());
	}


	/// Ritorna l'indice i di una variabile indicizzata del tipo "var(i). Se l'indice non c'e' torna -1"
	static int get_var_index(const string_t&	source)
	{
		// siccome questa funzione può essere chiamata molte volte, dichiaro regex di tipo statico
		// in modo che l'espressione venga compilata una sola volta per tutte
		static std::regex e("[\\s]*[[:alpha:]]+[a-zA-Z_0-9\\$\\.\\-]*\\([\\s]*([[:digit:]])+[\\s]*\\)[.]*");
		
		std::smatch result;
		string_t var_value;
		if (std::regex_search(source, result, e) && result.size()==2)
			return std::stoi(result[1]);
		else
			return -1;
	}


	/// inserisce una coppia var_name=var_val nella stringa dest. Il contenuto originale di dest viene cancellato
	///@ param dest la stringa di destinazione
	///@ var_name il nome della variabile
	///@ var_val  il valore della variabile
	template<typename T>
	static void put(string_t& dest, const string_t& var_name,  const T& var_val)
	{
		dest = var_name + equal() + to_string(var_val);
	}

		
	

	/**Aggiunge alla stringa dest la variabile var_name con il valore var_val 
	Se la variabile var_name e' gia' contenuta in dest, la sovrascrive
	param dest la stringa destinazione
	param var_name nome della variabile da aggiungere
	param var_val valore della variabile (int, double ,string,....)
	*/
	template<typename T>
	static void add(string_t& dest, const string_t& var_name,  const T& var_val)
	{
		// estraggo le variabili preesistenti dalla stringa di destinazione 
		var_set tmp_set=get(dest);
		// inserisco la nuova coppia nella mappa
		tmp_set[var_name]=to_string(var_val);
		// svuoto dest e poi inserisco il contenuto della mappa
		dest.clear();
		add(dest, tmp_set);
	}

	///Aggiunge alla stringa dest tutte le variabili(name, value) contenute  nella mappa 
	///associativa var_set
	///@param   dest      la stringa di destinazione
	///@param   var_set   una mappa associativa di coppie (var_name, var_value)
	static void add(string_t& dest, const var_set&  source_map)
	{
		var_set strstrmap=get(dest);
		strstrmap.insert(source_map.begin(), source_map.end());
		dest.clear();
		add(dest,strstrmap);
	}




	///Estrae il valore della variabile var_name dalla stringa source
	///@param   source     la stringa sorgente contenete le variabili
	///@param   var_name   il nome della variabile da estrarre
	///@param   var_value  il valore della variabile se questa esiste
	///@returns	 false se la variabile var_name non c'e' nella stringa source 
	///@returns	 true se la variabile e' stata estratta
	static bool get_var_value(const string_t& source, const string_t& var_name, string_t& var_value)
	{
		// 10-05-2011 nuovo
		// [\\s]* : 0 o +  spazi 
		std::smatch result;
		static std::regex e("[\\s]*[,]?\\<" + var_name + "[\\s]*=[\\s]*([a-zA-Z0-9_\\(\\)\\.\\-@]*)[\\s]*(?:[,;]|$)");
		if (std::regex_search(source, result, e))
		{
			var_value = result[1];
			return true;
		}
		else
			return false;
	}


	///Estrae tutte le variabili dalla stringa source e le copia nella mappa associativa var_set
	///@param    source      il record sorgente
	///@param    var_set   una mappa associativa di coppie (var_name, var_value)
	///@returns	 il numero di variabili trovate
	static int get(const string_t& source, var_set&  var_set)
	{
		std::regex e;
		//la prima espressione tra parentesi tonda contiene il nome della variabile
		//la seconda il valore in formato stringa con caratteri aggiuntivi { _ . - @} 
		//static STRING		name_format="[\\s]*([\\w\\.\\d]*)[\\s]*=[\\s]*([a-zA-Z0-9_ \\.\\-@\\s]*)";
		string_t		name_format=_T("[\\s]*([\\w\\.\\d]*)[\\s]*=[\\s]*([a-zA-Z0-9_\\.\\-@]*)[\\s]*");
		std::smatch	match;
		string_t		part;
		var_set.clear();
		e=name_format;
		for(string_t::size_type   p1=0, p2=0; p2!=string_t::npos; p1=p2+1)
		{
			p2=source.find(',',p1);
			// se non trova la virgola p2=npos per cui la substr(p1,p2-p1) tornerà la parte finale 
			part=source.substr(p1,p2-p1);
			if(std::regex_search(part,match, e) && match.size()==3)
					var_set[match[1]]=match[2];
		}
		return int(var_set.size());
	}




	///Estrae tutte le variabili dalla stringa source e le copia nella mappa associativa var_set che ritorna
	///@param       source      il record sorgente
	///@returns     una mappa associativa di coppie (var_name, var_value)
	static var_set get(const string_t& source)
	{
		var_set  v_set;
		get(source, v_set);
		return v_set;
	}
}; // end class






} // end namespaces PI
#endif