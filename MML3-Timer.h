#pragma once

#include<chrono>
#include<ctime>
#include<string>
#pragma warning(disable : 4996)
namespace MML3
{
	class Timer {
		typedef std::chrono::high_resolution_clock high_resolution_clock;
	public:
		explicit	Timer() :_start(std::chrono::high_resolution_clock::now()){}
		double 		start(){_start = std::chrono::high_resolution_clock::now();return now();}
		long long	nanoseconds() const{return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - _start).count();}
		double		seconds() const{return nanoseconds() / 1E9;}
		double		now() const{ return seconds(); }
		double		accuracy()const{	return 	double(std::chrono::high_resolution_clock::period::num)	/ std::chrono::high_resolution_clock::period::den ;	}

		static std::string current_time()
		{
			std::chrono::system_clock::time_point p = std::chrono::system_clock::now();
			std::time_t t = std::chrono::system_clock::to_time_t(p);
			return  std::ctime(&t);
		}
		
	private:
		std::chrono::high_resolution_clock::time_point _start;
	};







	
	


} // end namespace MML3



