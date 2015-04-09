#ifndef GNUPLOTTING_H
#define GNUPLOTTING_H

#include <string>
#include <cstdio>
#include <vector>
#include <armadillo>
using namespace std;

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

		void xrange(double &xmin, double &xmax);
		void yrange(double &ymin, double &ymax);
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
	private:
		string filename;
		FILE * pipe;
};

#endif

