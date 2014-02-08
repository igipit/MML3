

#include"PI-logstream.h"

namespace PI
{

	
	/*log_device_t log_device;
	log_stream_t log;
*/


	
	//std::ostream   log(& logbuf::I());
	std::ostream&   log		=	logbuf::ostream();
	std::wostream&  wlog = wlogbuf::ostream();


	
	


	std::string error(int id, const std::string& msg)	{ return "Error[" + std::to_string(id) + "](" + msg + "): "; }
	std::string warning(int id, const std::string& msg)	{
		return "Warning[" + std::to_string(id) + "](" + msg + ") : "; }
	std::string separator()	{ return "-----------------------------------------------------------\n"; }


}