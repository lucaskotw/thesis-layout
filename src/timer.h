#ifndef _TIMER_H
#define _TIMER_H


#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/date_time/microsec_time_clock.hpp>
#include <iostream>

class Timer
{
public:
    Timer(const std::string & name) : name(name),
        start(boost::date_time::microsec_clock<boost::posix_time::ptime>::local_time())
    {
    }

    ~Timer();

private:
    std::string name;
    boost::posix_time::ptime start;
};

#endif