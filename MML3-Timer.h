#pragma once

#include<chrono>
#include<ctime>
#include<string>
#include<vector>
#pragma warning(disable : 4996)

namespace MML3
{
    

struct Event
{
	Event(double t, const std::string& s) :time(t), name(s){}
	Event() = default;
	Event(const Event&) = default;
	double	   time = 0;
	std::string name; 
};
    
    class Timer {
	typedef std::chrono::high_resolution_clock              hr_clock;
	typedef std::chrono::high_resolution_clock::time_point	t_point;
	public:
	
		typedef std::vector<Event>			event_vector;

		explicit	Timer()			:start_(hr_clock::now()){ event_.reserve(10); new_event("start"); }
		double 		start()			{start_ = hr_clock::now(); event_.resize(0); return now(); }
		long long	nanoseconds()	{return std::chrono::duration_cast<std::chrono::nanoseconds>(hr_clock::now() - start_).count();}
		double		new_event(const std::string& name)
									{
										event_.push_back(Event(now(), name )); 
										return event_.back().time; 
									}

		
		double		seconds()		{return ((nanoseconds() / 1E9));}
		double		now()			{ return seconds(); }
		double		accuracy()const	{
										return 	double(std::chrono::high_resolution_clock::period::num)	/ std::chrono::high_resolution_clock::period::den ;
									}

		const event_vector& get_events()const{return event_;}

		/** Timer::current_time()
		* return a std::string with the current time
		*/
		static std::string current_time()
		{
			std::chrono::system_clock::time_point p = std::chrono::system_clock::now();
			std::time_t t = std::chrono::system_clock::to_time_t(p);
			return  std::ctime(&t);
		}
		
	private:
		t_point	start_;
		event_vector event_;
	};







	
	


} // end namespace MML3



