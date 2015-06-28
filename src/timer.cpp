#include "timer.h"


Timer::~Timer()
{
    using namespace std;
    using namespace boost;

    posix_time::ptime now(date_time::microsec_clock<posix_time::ptime>::local_time());
    posix_time::time_duration d = now - start;

    cout << name << " completed in " << d.total_milliseconds() / 1000.0 <<
        " seconds" << endl;
}