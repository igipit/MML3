
#include"PI-logstream.h"

namespace PI
{

	
	/*log_device_t log_device;
	log_stream_t log;
*/


	
	


	std::string error(int id, const std::string& msg)	{ return "Error[" + std::to_string(id) + "](" + msg + "): "; }
	std::string warning(int id, const std::string& msg)	{
		return "Warning[" + std::to_string(id) + "](" + msg + ") : "; }
	std::string separator()	{ return "-----------------------------------------------------------\n"; }

	std::string debug(const std::string& msg)	{ return "Debug [" + msg + "]\n"; }
	
}