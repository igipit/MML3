// log_stream.cpp : definisce il punto di ingresso dell'applicazione console.
//
#pragma once

#include<iomanip>
#include<string>
#include<vector>
#include <iostream> 


#define SOURCE_INFO "File:" << __FILE__ << "\tLine:" << __LINE__ << "\t"
//#define DEBUG_INFO "[DEBUG] File:" << __FILE__ << "\tLine:" << __LINE__ << "\t" << std::endl
#define PRINT_DEBUG_INFO std::cerr << "[DEBUG] File:" << __FILE__ << "\tLine:" << __LINE__ << "\t" << std::endl;
namespace PI
{


	std::string error(int id = 0, const std::string& msg ="");
	std::string warning(int id = 0, const std::string& msg="");
	std::string separator();
	std::string debug(const std::string& msg);


	template<typename CH>
	std::basic_ostream<CH>& get_stdio();

	template<>
	inline std::basic_ostream<char>& get_stdio<char>()
	{
		return std::cout;
	}

	template<>
	inline std::basic_ostream<wchar_t>& get_stdio<wchar_t>()
	{
		return std::wcout;
	}


	template<typename CH>
	class basic_logbuf : public std::basic_streambuf<CH>
	{
		typedef std::basic_ostream <CH> ostream_t;
	protected:
		
		// la definizione di questa funzione e' sufficente a far funzionare correttamente il buffer
		// in quanto le funzioni della classe basic_streambuf la utilizzeranno per l'output
		virtual int_type overflow(int_type c)
		{
			if (c != EOF)
			{
				bool eof = false;
				for (auto ps : os_)
				{
					ps->put(c);
					eof |= ps->eof();
				}
				if (eof)
					c = EOF;
			}
			return c;
		}
		// per efficienza sovrappongo anche la funzione che scrive piu' caratteri
		// se non viene definita, la versione di base utilizza overflow
		virtual std::streamsize xsputn(const CH* s, std::streamsize num)
		{

			bool failure = false;
			for (auto ps : os_)
			{
				ps->write(s,num);
				failure |= ps->fail();
			}
			if (failure)
				return 0;
			else
				return num;
		}

		// questa viene chiamata quando lo stream che possiede il buffer
		// richiede un flush
		virtual int sync()
		{
			bool failure = false;
			for (auto ps : os_)
			{
				ps->flush();
				failure |= ps->fail();
			}
			if (failure)
				return -1;
			else
				return 0;
		}

	public:
	
		~basic_logbuf()
		{
			// essendo gli oggetti a cui puntano le componenti di os_
			// gestiti esternamente non posso chiamare flush in quanto
			// quando viene chiamato il distruttore di basic_logbuf gli oggetti esterni 
			// potrebbero gia' essere distrutti
			for (auto ps : os_)
				ps=nullptr;
			
		}

	public:

		// static public interface

		/// Ritorna l'unica istanza di basic_logbuf
		static  basic_logbuf& I()
		{
			static basic_logbuf unique_instance(get_stdio<CH>());
			return unique_instance;
		}


		/// ritorna un output stream associato  all'unico  logbuf
		static ostream_t&  basic_logbuf::ostream()
		{
			static ostream_t os(&basic_logbuf::I());
			return os;
		}
		
		/// aggiunge uno stream di output al buffer
		static basic_logbuf& add_ostream(ostream_t& o)
		{

			I().os_.push_back(&o);
			return I();
		}

		/// rimuove un output stream dal buffer
		static void  remove_ostream(size_t i)
		{
			std::vector<ostream_t*>  os = I().os_;

			if (i < os.size())
			{
				os[i]->flush();
				os[i] = nullptr;;
				os.erase(os.begin() + i);
			}
		}

		/// ritorna il numero di output streams su cui il buffer riversa l'output
		static size_t num_ostreams(){ return I().os_.size(); }

	private:

		basic_logbuf(ostream_t& o){ os_.reserve(2); os_.push_back(&o); }
		std::vector<ostream_t*>  os_;

		

	};





	/// \brief This class is a derivate of basic_stringbuf which will output all the written data using the OutputDebugString function
	template<typename TChar, typename TTraits = std::char_traits<TChar>>
	class OutputDebugStringBuf : public std::basic_stringbuf<TChar, TTraits> {
	public:
		explicit OutputDebugStringBuf() : _buffer(256) {
			setg(nullptr, nullptr, nullptr);
			setp(_buffer.data(), _buffer.data(), _buffer.data() + _buffer.size());
		}

		~OutputDebugStringBuf() {
		}

		static_assert(std::is_same<TChar, char>::value || std::is_same<TChar, wchar_t>::value, "OutputDebugStringBuf only supports char and wchar_t types");

		int sync() try {
			MessageOutputer<TChar, TTraits>()(pbase(), pptr());
			setp(_buffer.data(), _buffer.data(), _buffer.data() + _buffer.size());
			return 0;
		}
		catch (...) {
			return -1;
		}

		int_type overflow(int_type c = TTraits::eof()) {
			auto syncRet = sync();
			if (c != TTraits::eof()) {
				_buffer[0] = c;
				setp(_buffer.data(), _buffer.data() + 1, _buffer.data() + _buffer.size());
			}
			return syncRet == -1 ? TTraits::eof() : 0;
		}


	private:
		std::vector<TChar>		_buffer;

	};







	
	typedef basic_logbuf<char>		logbuf;
	typedef basic_logbuf<wchar_t>	wlogbuf;
	
	// log and wlog are provided for user convenience since it is more beatiful to write
	// log << "log message";
	//    instead of
	// logbuf::ostream() << "log message";
	//
	// but they are really global references to the unique objects
	// basic_logbuf<char>::ostream()
	// and, respectively, to
	// basic_logbuf<wchar_t>::ostream()
	extern std::ostream&			log;
	extern std::wostream&			wlog;
	



	/** Usage
	when you include PI-logstream.cpp in the project you have a global object (one and only one) PI::log (and PI::wlog)
	of type std::ostream& (and std::wostream&) for sending the log output (by default it redirects to the  stdio terminal)

	if you need a multiple output add  new ostreams to the log buffer by means of static methods
	PI::logbuf::add(new_ostream1);
	PI::logbuf::add(new_ostream2);
	...
	if the added ostream is connected to a file (if it is an ofstream) it is your responsibility 
	to open and close the stream
	
	*/














} // end namespace Pi



