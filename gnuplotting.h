#ifndef GNUPLOTTING_H
#define GNUPLOTTING_H

#include <string>
#include <cstdio>
#include <vector>
#include <armadillo>
using namespace std;
using namespace arma;

/*
Gnuplotting & operator << (Gnuplotting & outstream, char * command){
	outstream.cmd(command);
	return outstream;
}
*/

class Gnuplotting {
	public:
		Gnuplotting();
		Gnuplotting(string filename);
		Gnuplotting(char *filename);
		~Gnuplotting();

		void xrange(double xmin, double xmax);
		void yrange(double ymin, double ymax);
		void title(string & titlename);
		void title(const char * titlename);
		void xlabel(string & x);
		void ylabel(string & y);
		void xlabel(const char * x);
		void ylabel(const char * y);
		void cmd(string & command);
		void cmd(const char *  command);

		Gnuplotting & operator<<(const char * command);

		void xystream(vector<double> & x, vector<double> & y);
		void xystream(size_t & N, arma::vec & x, arma::vec & y);
		void xystream_replot(size_t & N, arma::vec & x, arma::vec & y, const char* title);
		void two_xystream(size_t &N1, vec &x1, vec &y1, const char* title1, size_t &N2, vec &x2, vec &y2, const char* title2);
		void xyzstream(size_t & N, arma::vec & x, arma::vec & y, arma::mat & z);
		void heatmap(size_t & Nx, size_t & Ny, arma::mat & z);
	private:
		string filename;
		FILE * pipe;
};

#endif

