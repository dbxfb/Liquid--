#ifndef TIMER_HH_
#define TIMER_HH_

class Timer
{
public:
    static u64 getTimeUs();

public:
    Timer();
    ~Timer();

    f64 delta();

private:
    u64 mLastTime;
};


#endif /* TIMER_HH_ */
