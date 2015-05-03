#include "diffusion1d.h"
#include <cmath>
#include <gsl/gsl_sf_erf.h>
#include <iostream>
using namespace arma;
using namespace std;

Diffusion1d::Diffusion1d(){

theta = 0.5;
u_0 = 1.0;
L = 1.0;
D_0 = 1.0;
tau = L*L/D_0;
isAbsorbing = true;
isReflective = false;
t = 0.0;
N = 101;
dt = 0.0001;
dx = 1.0/((double) (N-1));
x = zeros<vec>(N);
u = zeros<vec>(N);
u_unbounded = zeros<vec>(N);
u_absorbing = zeros<vec>(N);
u_reflective = zeros<vec>(N);

createX();
createDeltaFunction();
}

Diffusion1d::~Diffusion1d(){

}

void Diffusion1d::createDeltaFunction(){
	u = zeros<vec>(N);
	bool even = false;
	if(N % 2 == 0){
		even = true;
	}
	if(even){
		size_t right = N/2;
		size_t left = N/2 - 1;
		double height = u_0/(2.0*dx);
		u(right) = height;
		u(left) = height;
		u_peak = height;
	}else{
		size_t mid = (N+1)/2 - 1;
		size_t right = mid + 1;
		size_t left = mid - 1;
		double height = u_0/(3.0*dx);
		u(mid) = height;
		u(right) = height;
		u(left) = height;
		u_peak = height;
	}
	double ymax = u_peak;
}

void Diffusion1d::createCustomInitial(){
	for(size_t i = 0; i<N; i++){
		u(i) = u_initial(x(i));
	}
}

void Diffusion1d::setReflective(){
	isReflective = true;
	isAbsorbing = false;
}

double Diffusion1d::getTime(){
	return t;
}

void Diffusion1d::setAbsorbing(){
	isReflective = false;
	isAbsorbing = true;
}

void Diffusion1d::setL(double & l){ L = l; }

void Diffusion1d::setD_0(double & d_0){ D_0 = d_0; }

double Diffusion1d::alpha(size_t i){
	return f(x(i))*dt/(dx*dx);
}

double Diffusion1d::beta(size_t i){
	return f_prime(x(i))*dt/(2.0*dx);
}


void Diffusion1d::createDiagonals(){
	if(isReflective){
		B_diag_0 = zeros<vec>(N);
		B_diag_up = zeros<vec>(N-1);
		B_diag_down = zeros<vec>(N-1);
		for(size_t i = 0; i < N; i++){
			B_diag_0(i) = 1.0 - 2.0*(1.0-theta)*alpha(i);
		}
		B_diag_up(0) = 2.0*(1.0-theta)*alpha(0);
		for(size_t i = 1; i < N-1; i++){
			B_diag_up(i) = (1.0-theta)*(alpha(i) + beta(i));
		}
		for(size_t i = 0; i < N-2; i++){
			B_diag_down(i) = (1.0-theta)*(alpha(i+1) - beta(i+1));
		}
		B_diag_down(N-2) = 2.0*(1.0-theta)*alpha(N-1);
	}
	else{
		B_diag_0 = zeros<vec>(N-2);
		B_diag_up = zeros<vec>(N-3);
		B_diag_down = zeros<vec>(N-3);
		for(size_t i = 0; i < N-2; i++){
			B_diag_0(i) = 1.0 - 2.0*(1.0-theta)*alpha(i+1);
		}
		for(size_t i = 0; i < N-3; i++){
			B_diag_up(i) = (1.0-theta)*(alpha(i+1) + beta(i+1));
			B_diag_down(i) = (1.0-theta)*(alpha(i+2) - beta(i+2));
		}
	}
}

void Diffusion1d::createB(){
	if(isReflective){
		B = zeros<mat>(N,N);
	}
	else{
		B = zeros<mat>(N-2,N-2);
	}
	B.diag(0) = B_diag_0;
	B.diag(-1) = B_diag_down;
	B.diag(1) = B_diag_up;
}

void Diffusion1d::createX(){
	for(size_t i = 0; i < N; i++){
		x(i) = ((double) i) * dx;
	}
}

void Diffusion1d::createUnbounded_analytic(){
	for(size_t i = 0; i<N; i++){
		u_unbounded(i) = (u_0/sqrt(4*M_PI*t)*exp(-((x(i)-0.5)*(x(i)-0.5))/(4.0*t)));
	}
}

double Diffusion1d::u_n_eigen_absorbing(double x, size_t n){
	if(n==0){
		return 0.0;
	}else{
		double n_double = (double) n;
		return sqrt(2.0)*sin(n_double*M_PI*x);
	}
}

double Diffusion1d::u_n_eigen_reflective(double x, size_t n){
	if(n==0){
		return 1.0;
	}else{
		double n_double = (double) n;
		return sqrt(2.0)*cos(n_double*M_PI*x);
	}
}

void Diffusion1d::createAbsorbing_analytic(size_t m){
	u_absorbing = zeros<vec>(N);
	for(size_t i = 0; i<N; i++){
		for(size_t n = 0; n<m; n++){
			u_absorbing(i) = u_absorbing(i) + exp(-(n*n*M_PI*M_PI*t))*u_n_eigen_absorbing(x(i),n)*u_n_eigen_absorbing(0.5,n);
		}
			}
	u_absorbing = u_0*u_absorbing;
}

void Diffusion1d::createReflective_analytic(size_t m){
	u_reflective = zeros<vec>(N);
	for(size_t i = 0; i<N; i++){
		for(size_t n = 0; n<m; n++){
			u_reflective(i) = u_reflective(i) + exp(-(n*n*M_PI*M_PI*t))*u_n_eigen_reflective(x(i),n)*u_n_eigen_reflective(0.5,n);
		}
			}
	u_reflective = u_0*u_reflective;
}

//PLOTTING
void Diffusion1d::plot(){
	gplt.xystream(N,x,u);
}

void Diffusion1d::plot_with_unbounded(){
	gplt.two_xystream(N,x,u,"Numerical",N,x,u_unbounded,"Unbounded analytical");
}

void Diffusion1d::plot_with_absorbing(){
	gplt.two_xystream(N,x,u,"Numerical",N,x,u_absorbing,"Absorbing analytical");
}

void Diffusion1d::plot_with_reflective(){
	gplt.two_xystream(N,x,u,"Numerical",N,x,u_reflective,"Reflective analytical");
}

void Diffusion1d::plotUnbounded(){
	gplt.xystream(N,x,u_unbounded);
}

void Diffusion1d::plotAbsorbing(){
	gplt.xystream(N,x,u_absorbing);
}

void Diffusion1d::plotReflective(){
	gplt.xystream(N,x,u_reflective);
}

void Diffusion1d::initialize(){
	t = 0.0;
	createDiagonals();
	createB();
	createDeltaFunction();
}

void Diffusion1d::solve(){
	if(isReflective){
		mat A = zeros<mat>(N,N);
		vec A_diag_0 = zeros<vec>(N);
		vec A_diag_up = zeros<vec>(N-1);
		vec A_diag_down = zeros<vec>(N-1);
		for(size_t i = 0; i < N; i++){
			A_diag_0(i) = 1.0 + 2.0*theta*alpha(i);
		}
		A_diag_up(0) = -2.0*theta*alpha(0);
		for(size_t i = 1; i < N-1; i++){
			A_diag_up(i) = -theta*(alpha(i) + beta(i));
		}
		for(size_t i = 0; i < N-2; i++){
			A_diag_down(i) = theta*(beta(i+1) - alpha(i+1)) ;
		}
		A_diag_down(N-2) = -2.0*theta*alpha(N-1);
		A.diag(0) = A_diag_0;
		A.diag(-1) = A_diag_down;
		A.diag(1) = A_diag_up;
		A_inv = inv(A);
	}
	else{
		mat A = zeros<mat>(N-2,N-2);
		vec A_diag_0 = zeros<vec>(N-2);
		vec A_diag_up = zeros<vec>(N-3);
		vec A_diag_down = zeros<vec>(N-3);
		for(size_t i = 0; i < N-2; i++){
			A_diag_0(i) = 1.0 + 2.0*theta*alpha(i+1);
		}
		for(size_t i = 0; i < N-3; i++){
			A_diag_up(i) = -theta*(alpha(i+1) + beta(i+1));
			A_diag_down(i) = theta*(beta(i+2) - alpha(i+2));
		}
		A.diag(0) = A_diag_0;
		A.diag(-1) = A_diag_down;
		A.diag(1) = A_diag_up;
		A_inv = inv(A);
	}
}

void Diffusion1d::iterate(size_t it){
	if(isReflective){
		for(size_t i = 0; i<it; i++){
			u = A_inv*B*u;
			t = t + dt;
		}
	}else{
		vec temp = u.subvec(1,N-2);
		for(size_t i = 0; i<it; i++){
			temp = A_inv*B*temp;
			t = t + dt;
		}
		u.subvec(1,N-2) = temp;
	}
}

ConstantDiffusion::ConstantDiffusion(): Diffusion1d(){
	initialize();
}


double ConstantDiffusion::f(double & x){
	return 1.0;
}

double ConstantDiffusion::f_prime(double & x){
	return 0.0;
}

double ConstantDiffusion::u_initial(double & x){
	return 1.0;
}


StepDiffusion::StepDiffusion(): Diffusion1d(){
//	initialize();
	gamma_plus = 1.0;
	gamma_minus = 1.0;
	u_step = zeros<vec>(N);
	initialize();

}


double StepDiffusion::f(double & x){
	if(x < 0.5){
		return gamma_minus;
	}else{
		return gamma_plus;
	}
}

double StepDiffusion::f_prime(double & x){
	return 0.0;
}

double StepDiffusion::u_initial(double & x){
	return 1.0;
}

void StepDiffusion::createStep_analytic(){
	u_step = zeros<vec>(N);
	//double A_plus = 2.0/(1.0+gsl_sf_erf(0.5/(sqrt(4.0*gamma_plus*t))) + sqrt(gamma_minus/gamma_plus)*
	//		exp((gamma_plus-gamma_minus)*0.5*0.5/(4.0*gamma_plus*gamma_minus*t))*(1.0-gsl_sf_erf(0.5/(sqrt(4.0*gamma_minus*t)))));
	double A_plus = 2.0/(1.0 + sqrt(gamma_minus/gamma_plus));
	cout << t;
	vec y = x - 0.5;
	for(size_t i=0; i<N; i++){
		if(y(i) < 0){
			//double A_minus = sqrt(gamma_minus/gamma_plus)*exp((gamma_plus-gamma_minus)*0.5*0.5/(4.0*gamma_plus*gamma_minus*t))*A_plus;
			double A_minus = A_plus*sqrt(gamma_minus/gamma_plus);
			u_step(i) = (A_minus/(sqrt(4.0*M_PI*gamma_minus)))*exp(-(y(i)*y(i))/(4.0*gamma_minus*t));
		}else{
			u_step(i) = (A_plus/(sqrt(4.0*M_PI*gamma_plus)))*exp(-(y(i)*y(i))/(4.0*gamma_plus*t));
		}
	}
	u_step=u_step*10.0;
}

void StepDiffusion::plot_with_step(){
	gplt.two_xystream(N,x,u,"Numerical",N,x,u_step,"Step analytical");
}

void StepDiffusion::plot_step(){
	gplt.xystream(N,x,u_step);
}

LinearDiffusion::LinearDiffusion(): Diffusion1d(){
	initialize();
}

double LinearDiffusion::f(double & x){
	return 1.0;
}

double LinearDiffusion::f_prime(double & x){
	return 0.0;
}

double LinearDiffusion::u_initial(double & x){
	return 1.0;
}

SinusDiffusion::SinusDiffusion(): Diffusion1d(){
	initialize();
}

double SinusDiffusion::f(double & x){
	return sin(M_PI*x);
}

double SinusDiffusion::f_prime(double & x){
	return M_PI*cos(M_PI*x);
}

double SinusDiffusion::u_initial(double & x){
	return 1.0;
}

SawDiffusion::SawDiffusion(): Diffusion1d(){
	initialize();
}

double SawDiffusion::f(double & x){
	return 1;
}

double SawDiffusion::f_prime(double & x){
	return 0.0;
}

double SawDiffusion::u_initial(double & x){
	double p = 5.0;
	double l = 1.0/p;
	double y = x/l;
	double decimal;
	double intpart;
	decimal = modf(y,&intpart);
	double z = decimal - intpart*l;
	if(decimal < 0.5 ){
		return 2.0*z/l;
	}else{
		return (2.0-2.0*z/l);
	}
}
