#include "diffusion1d.h"
#include <unistd.h>
#include <iostream>

using namespace std;
using namespace arma;

int main(){
/*
ConstantDiffusion absorb, reflect;
absorb.solve();
reflect.setReflective();
reflect.initialize();
reflect.solve();

//absorb.iterate(1);
//reflect.iterate(1);

for(size_t i=0; i < 10; i++){
	absorb.iterate(1);
	absorb.createAbsorbing_analytic(1000);
	reflect.iterate(1);
	reflect.createReflective_analytic(1000);
	
	absorb.plot_with_absorbing();
	reflect.plot_with_reflective();
	sleep(1);
}
*/

StepDiffusion step;
step.setReflective();
step.initialize();
step.solve();

for(size_t i=0; i < 1200; i++){
	step.iterate(4);
	step.plot();
	usleep(1e4);
}

LinearDiffusion linear;
linear.setReflective();
linear.initialize();
linear.createCustomInitial();
linear.solve();

for(size_t i=0; i < 200; i++){
	linear.iterate(2);
	linear.plot();
	usleep(1e4);
}

SinusDiffusion sinus;
sinus.setReflective();
sinus.initialize();
sinus.solve();

for(size_t i=0; i < 200; i++){
	sinus.iterate(1);
	sinus.plot();
	usleep(1e4);
}

/*
absorb.createReflective_analytic(1000);
absorb.createAbsorbing_analytic(1000);
absorb.plotAbsorbing();
sleep(3);
absorb.plotReflective();
sleep(10);
*/
/*
for(size_t i=0; i < 200; i++){
	absorb.iterate(5);
	absorb.createUnbounded_analytic();
	absorb.plot_with_unbounded();
	usleep(1e4);
}
*/
/*
for(size_t i = 0; i < 300; i++){
	absorb.iterate(10);
	reflect.iterate(10);
	absorb.createUnbounded_analytic();
	reflect.createUnbounded_analytic();
	absorb.plot_with_absorbing();
	reflect.plot_with_reflective();

	usleep(1e5);
}
*/
	return 0;
}
