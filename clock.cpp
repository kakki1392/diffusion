#include "clock.h"
using namespace std;

const uint64_t billion = 1000000000;
const uint64_t million = 1000000;
const uint64_t thousand = 1000;

Clock::Clock(){
	clock_gettime(CLOCK_MONOTONIC, &time);
}

Clock::~Clock(){
}

void Clock::start(){
	clock_gettime(CLOCK_MONOTONIC, &time);
}

time_t Clock::stop(){
	timespec temp;
	clock_gettime(CLOCK_MONOTONIC, &temp);
	return temp.tv_sec - time.tv_sec;
}

uint64_t Clock::mstop(){
	timespec temp;
	clock_gettime(CLOCK_MONOTONIC, &temp);
	uint64_t ndiff = temp.tv_nsec - time.tv_nsec;
	return (billion*(temp.tv_sec - time.tv_sec) + ndiff)/million;
}

uint64_t Clock::ustop(){
	timespec temp;
	clock_gettime(CLOCK_MONOTONIC, &temp);
	uint64_t ndiff = temp.tv_nsec - time.tv_nsec;
	return (billion*(temp.tv_sec - time.tv_sec) + ndiff)/thousand;
}

uint64_t Clock::nstop(){
	timespec temp;
	clock_gettime(CLOCK_MONOTONIC, &temp);
	uint64_t ndiff = temp.tv_nsec - time.tv_nsec;
	return (billion*(temp.tv_sec - time.tv_sec) + ndiff);
}

void Clock::wait(time_t sec){
	timespec start;
	timespec stop;
	clock_gettime(CLOCK_MONOTONIC, &start);
	do{
		clock_gettime(CLOCK_MONOTONIC, &stop);
	}while((stop.tv_sec - start.tv_sec) < sec);
}
	
void Clock::mwait(uint64_t msec){
	timespec start;
	timespec stop;
	uint64_t ndiff;
	uint64_t mdiff;
	clock_gettime(CLOCK_MONOTONIC, &start);
	do{
		clock_gettime(CLOCK_MONOTONIC, &stop);
		ndiff = stop.tv_nsec - start.tv_nsec;
		mdiff = (billion*(stop.tv_sec - stop.tv_sec) + ndiff)/million;
	}while(mdiff < msec);
}

void Clock::uwait(uint64_t usec){
	timespec start;
	timespec stop;
	uint64_t ndiff;
	uint64_t udiff;
	clock_gettime(CLOCK_MONOTONIC, &start);
	do{
		clock_gettime(CLOCK_MONOTONIC, &stop);
		ndiff = stop.tv_nsec - start.tv_nsec;
		udiff = (billion*(stop.tv_sec - stop.tv_sec) + ndiff)/thousand;
	}while(udiff < usec);
}

void Clock::nwait(uint64_t nsec){
	timespec start;
	timespec stop;
	uint64_t ndiff;
	uint64_t nndiff;
	clock_gettime(CLOCK_MONOTONIC, &start);
	do{
		clock_gettime(CLOCK_MONOTONIC, &stop);
		ndiff = stop.tv_nsec - start.tv_nsec;
		nndiff = billion*(stop.tv_sec - stop.tv_sec) + ndiff;
	}while(nndiff < nsec);
}


