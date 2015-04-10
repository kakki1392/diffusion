#ifndef CLOCK_H
#define CLOCK_H
#include <time.h>
#include <stdint.h>
using namespace std;

/*
  Class for high resolution timing.
  start(): FIRST LOOK AT CLOCK
  *stop(): NEXT LOOK AT CLOCK, THIS MANY *SECONDS HAVE PASSED
  *wait(*SECONDS): WAIT FOR *SECONDS, THIS WILL PAUSE PROGRAM

  "*" is one of: " ", "m", "u", "n".

  " " : none
  "m" : milli
  "u" : micro
  "n" : nano
 */
class Clock{
	public:
		Clock();
		~Clock();
		void start();
		time_t stop();
		uint64_t mstop();
		uint64_t ustop();
		uint64_t nstop();
		void wait(time_t sec);
		void mwait(uint64_t msec);
		void uwait(uint64_t usec);
		void nwait(uint64_t nsec);
	private:
		timespec time;
};

#endif
