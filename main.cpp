#include "diffusion1d.h"
#include <unistd.h>
#include <iostream>

using namespace std;
using namespace arma;

int main(){

ConstantDiffusion A;
A.solve();

for(size_t i = 0; i < 200; i++){
	A.iterate(2);
	A.plot();
	usleep(1e5);
}
/*
A.setReflective();
A.initialize();
A.solve();
for(size_t i = 0; i < 200; i++){
	A.iterate(2);
	A.plot();
	usleep(1e5);
}
*/
	return 0;
}
