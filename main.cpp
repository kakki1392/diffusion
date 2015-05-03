#include "diffusion1d.h"
#include <unistd.h>
#include <iostream>
#include <sstream>

using namespace std;
using namespace arma;

int main(){
	ConstantDiffusion A, B;
	A.solve();
	B.setReflective();
	B.initialize();
	B.solve();

	A.iterate(20);
	B.iterate(20);
	double t1 = A.getTime();
	double t2 = B.getTime();
	stringstream s1, s2;
	s1 << "set title 'Diffusion, absorbing boundary, t = " << t1 << "'";
	s2 << "set title 'Diffusion, reflective boundary, t = " << t2 << "'";
	string title1 = s1.str();
	string title2 = s2.str();

	A.createAbsorbing_analytic(100);
	B.createReflective_analytic(100);


 	A.gplt.cmd("set term pdfcairo");
	A.gplt.cmd("set output 'absorbing1.pdf'");
	B.gplt.cmd("set term pdfcairo");
	B.gplt.cmd("set output 'reflective1.pdf'");

	//A.gplt.cmd("set title 'test1'");
	A.gplt.cmd(title1);
	A.gplt.cmd("set xlabel 'x'");
	A.gplt.cmd("set ylabel 'u(x,t)'");
	A.gplt.cmd("set yrange [0:7]");

	
	//B.gplt.cmd("set title 'test1'");
	B.gplt.cmd(title2);
	B.gplt.cmd("set xlabel 'x'");
	B.gplt.cmd("set ylabel 'u(x,t)'");
	B.gplt.cmd("set yrange [0:7]");

//	A.gplt.cmd("set multiplot layout 2,2");
//	B.gplt.cmd("set multiplot layout 2,2");

	A.plot_with_absorbing();
	B.plot_with_reflective();

	A.gplt.cmd("unset output");
	B.gplt.cmd("unset output");


	StepDiffusion C;
	C.solve();
	C.iterate(10);
	C.createStep_analytic();
	C.plot_with_step();
/*
	for(int i=2; i < 5; i++){
		A.iterate(100);
		B.iterate(100);
		A.createAbsorbing_analytic(100);
		B.createReflective_analytic(100);
		
		A.plot_with_absorbing();
		B.plot_with_reflective();
		sleep(1);
	}
	//A.gplt.cmd("unset multiplot");
	//B.gplt.cmd("unset multiplot");

	A.gplt.cmd("unset output");
	B.gplt.cmd("unset output");

*/

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
/*
StepDiffusion step;
step.setReflective();
step.initialize();
step.solve();
//step.gplt.cmd("set yrange [0:100]");
for(size_t i=0; i<100; i++){
	step.iterate(1);
	step.plot();
	usleep(1e5);
}
*/
/*
for(size_t i=0; i < 1500; i++){
	step.iterate(30);
	step.plot();
	usleep(1e4);
}
*/
/*
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
*/

/*
SinusDiffusion sinus;
sinus.setReflective();
sinus.initialize();
sinus.solve();

for(size_t i=0; i < 300; i++){
	sinus.iterate(4);
	sinus.plot();
	usleep(1e4);
}
*/
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
