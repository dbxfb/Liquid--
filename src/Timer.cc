#include "common.hh"

#include <sys/time.h>

/* ************************************************************************** *
 * Timer : Static methods                                                     *
 * ************************************************************************** */
u64 Timer::getTimeUs()
{
    timeval tv;
    gettimeofday(&tv, NULL);

    return u64(tv.tv_sec) * 1000000ULL + u64(tv.tv_usec);
}

/* ************************************************************************** *
 * Timer : Public methods                                                     *
 * ************************************************************************** */

Timer::Timer()
{
    mLastTime = getTimeUs();
}

Timer::~Timer()
{

}

f64 Timer::delta()
{
    u64 time = getTimeUs();
    f64 d = f64(time - mLastTime) / 1000000.0;
    mLastTime = time;

    return d;
}
