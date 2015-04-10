#include "clock.h"
#include <iostream>
#include <unistd.h>
#include <time.h>

using namespace std;

int main(){
	Clock clock;
	clock.start();
	sleep(1);
	cout << clock.stop() << endl;
	cout << clock.mstop() << endl;
	cout << clock.ustop() << endl;
	cout << clock.nstop() << endl << endl;

	clock.start();
	clock.wait(3);
	cout << clock.stop() << endl;
	cout << clock.mstop() << endl;
	cout << clock.ustop() << endl;
	cout << clock.nstop() << endl << endl;

	clock.start();
	sleep(2);
	cout << clock.nstop() << endl;




	return 0;
}
