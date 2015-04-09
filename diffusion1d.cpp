#include "diffusion1d.h"
using namespace arma;
using namespace std;


Diffusion1d::Diffusion1d(){

theta = 0.5;
u_0 = 100.0;
L = 1.0;
D_0 = 1.0;
tau = L*L/D_0;
isAbsorbing = false;
isReflective = true;
t = 0.0;
N = 70;
dt = 0.0001;
dx = 1.0/((double) (N-1));
x = zeros<vec>(N);
u = zeros<vec>(N);
createX();
createDeltaFunction();
}

Diffusion1d::~Diffusion1d(){

}

void Diffusion1d::createDeltaFunction(){
	bool even = false;
	if(N % 2 == 0){
		even = true;
	}
	if(even){
		size_t right = N/2;
		size_t left = N/2 - 1;
		u(right) = u_0/dx;
		u(left) = u_0/dx;
	}else{
		size_t mid = (N+1)/2 - 1;
		size_t right = mid + 1;
		size_t left = mid - 1;
		u(mid) = u_0/(2.0*dx);
		u(right) = u_0/(2.0*dx);
		u(left) = u_0/(2.0*dx);
	}
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

void Diffusion1d::setAbsorbing(){
	isReflective = false;
	isAbsorbing = true;
}

void Diffusion1d::setL(double & l){ L = l; }

void Diffusion1d::setD_0(double & d_0){ D_0 = d_0; }

void Diffusion1d::print_dt_SI(){
	cout << dt*tau << endl;
}

void Diffusion1d::print_dx_SI(){
	cout << dx*L << endl;
}

void Diffusion1d::print_x(){
	x.print("x: ");
}

void Diffusion1d::print_u(){
	u.print("u: ");
}

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

void Diffusion1d::initialize(){
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

void Diffusion1d::print_AB(){
	mat temp = A_inv*B;
	temp.print();
}

void Diffusion1d::print_B(){
	B.print("B: ");
}

void Diffusion1d::print_A_inv(){
	A_inv.print("A_inv: ");
}

void Diffusion1d::plot(){
	double ymin = 0.0;
	double ymax = 1500.0;
	gplt.yrange(ymin,ymax);
	gplt.xystream(N,x,u);
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





