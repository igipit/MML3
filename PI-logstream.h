// log_stream.cpp : definisce il punto di ingresso dell'applicazione console.
//
#pragma once

#include<iomanip>
#include<string>
#include <iostream> 
#include<fstream>


//#define SOURCE_INFO0 "File:"
//#define SOURCE_INFO1 SOURCE_INFO0 ## __FILE__
//#define SOURCE_INFO2 SOURCE_INFO1 ##"\tLine:"
//#define SOURCE_INFO3 SOURCE_INFO2 ##__LINE__
//#define SOURCE_INFO  SOURCE_INFO3 ##"\t"

//#define DEBUG_INFO "[DEBUG] File:" << __FILE__ << "\tLine:" << __LINE__ << "\t" << std::endl
#define PRINT_DEBUG_INFO std::cerr << "[DEBUG] File:" << __FILE__ << "\tLine:" << __LINE__ << "\t" << std::endl;
namespace PI
{


	std::string error(int id = 0, const std::string& msg ="");
	std::string warning(int id = 0, const std::string& msg="");
	std::string separator();
	std::string debug(const std::string& msg);

	

	template<typename CH>
	class basic_logbuf : public std::basic_streambuf<CH>
	{
		typedef std::basic_ostream <CH>		ostream_t;
		typedef std::basic_ofstream<char>	ofstream_t;
	public:
		
		// la definizione di questa funzione e' sufficente a far funzionare correttamente il buffer
		// in quanto le funzioni della classe basic_streambuf la utilizzeranno per l'output
		virtual int_type overflow(int_type c)
		{
			if (c != EOF)
			{
				bool eof = false;
				if (os_)
				{
					os_->put(c);
					eof|= os_->eof();
				}
				if (file_stream_.is_open())
				{
					file_stream_.put(c);
					eof |= file_stream_.eof();
				}
				
				if (eof)
					return  EOF;
			}
			return c;
		}
		// per efficienza sovrappongo anche la funzione che scrive piu' caratteri
		// se non viene definita, la versione di base utilizza overflow
		virtual std::streamsize xsputn(const CH* s, std::streamsize num)
		{

			bool failure = false;
			if (os_)
			{
				os_->write(s,num);
				failure |= os_->fail();
			}

			if (file_stream_.is_open())
			{
				file_stream_.write(s, num);
				failure |= file_stream_.fail();
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
			
			if(os_)
			{
				os_->flush();
				failure |= os_->fail();
			}
			if(file_stream_.is_open())
			{
				file_stream_.flush();
				failure |= file_stream_.fail();
			}
			if (failure)
				return -1;
			else
				return 0;
		}

	
	
		~basic_logbuf()
		{
			if (file_stream_.is_open())
				file_stream_.close();
		}


		
		bool open(const std::string& filename, std::ios_base::openmode mode = std::ios_base::out)
		{
			close();
			file_stream_.open(filename, mode);
			return file_stream_.is_open();
		}

		void close()
		{
			if (file_stream_.is_open())
				file_stream_.close();
		}

		basic_logbuf(ostream_t& o ):os_(&o){}

		ostream_t *  os_ =nullptr;
		ofstream_t	file_stream_;
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
	



	/** Usage
	when you include PI-logstream.cpp in the project you have a global object (one and only one) std::clog (and PI::wlog)
	of type std::ostream& (and std::wostream&) for sending the log output (by default it redirects to the  stdio terminal)

	if you need a multiple output add  new ostreams to the log buffer by means of static methods
	std::clogbuf::add(new_ostream1);
	std::clogbuf::add(new_ostream2);
	...
	if the added ostream is connected to a file (if it is an ofstream) it is your responsibility 
	to open and close the stream
	
	*/














} // end namespace Pi



