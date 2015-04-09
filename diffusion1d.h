#ifndef DIFFUSION1D_H
#define DIFFUSION1D_H
#include <armadillo>
#include "gnuplotting.h"

using namespace arma;
using namespace std;

/* One dimensional diffusion in 1D. This is an abstract class, only specialized systems 
 * may be instanced. Diffusion constant may vary with position x. Either
 * reflective or absorbing boundaries. Initial distribution is calculated by u_initial(x).
 * Inherited classes will have own implementations of initial
 * distribution and the variation of the diffusion constant.
 */

/* Diffusion problem is del_t u = del_x D(x) del_x u.
 * u -> u_vector, where u_vector contains all u=u(x), x is discretized.
 * Typical numerical implementation: A_matrix * u_vector_(n+1) = B_matrix * u_vector_(n),
 * where (n+1) is the next time step after (n).
 * A_matrix must be inverted, and then the system can be iterated through time.
 * The numerical scheme is the theta-scheme of order O(dx^2) O(dt^2), generalized
 * Crank-Nicholson.
 */

class Diffusion1d{
	public:
		Diffusion1d();
		~Diffusion1d();
		virtual double f(double & x) = 0;
		virtual double f_prime(double & x) = 0;
		virtual double u_initial(double & x) = 0;

		void createCustomInitial();

		void setReflective();
		void setAbsorbing();
		void setL(double & l);
		void setD_0(double & d_0);

		void plot();
		
		void print_dt_SI();
		void print_dx_SI();
		void print_B();
		void print_A_inv();
		void print_AB();
		void print_x();
		void print_u();

		void initialize();    //Initializes matrices and vectors 
		void solve();         //Computes A_matrix_inverse
		void iterate(size_t it);       //Iterates the system through time
		void iterateForSeconds(double seconds);

	private:
		double alpha(size_t i);
		double beta(size_t i);
		void createDiagonals();
		void createB();
		void createX();
		void createDeltaFunction();

		double u_0; //number of particles
		double L;   //size of system
		double D_0; //max size of diffusion constant D=D(x)
		double tau; //characteristic time scale, tau=D_0/L*L

		Gnuplotting gplt;

		bool isReflective;
		bool isAbsorbing;

		double theta;
		double t;  //current time of system
		size_t N;  //Number of spatial points in [0,1]
		double dx; //x=jdx, j=0,....N-1
		double dt; //

		vec x;  //arma::vector containing all positions
		vec u;  //arma::vector containing u(x) at time t

		mat A_inv;  //inverse of A_matrix, is only inverted ONCE
		mat B;
		vec B_diag_0;  //diagonals of B_matrix
		vec B_diag_up;
		vec B_diag_down;
};

class ConstantDiffusion: public Diffusion1d{
	public:
		ConstantDiffusion();
		~ConstantDiffusion(){};
		virtual double f(double & x);
		virtual double f_prime(double & x);
		virtual double u_initial(double & x);

};

#endif
