#include "gnuplotting.h"
#include <cstdio>
#include <string>
#include <sstream>
#include <armadillo>

using namespace std;

Gnuplotting::Gnuplotting(){
	pipe = popen("gnuplot -persist", "w");
}

Gnuplotting::~Gnuplotting(){
	pclose(pipe);
}

void Gnuplotting::cmd(const char * command){
	fprintf(pipe,"%s\n",command);
	fflush(pipe);
}

void Gnuplotting::cmd(string & command){
	fprintf(pipe,"%s\n",command.c_str());
	fflush(pipe);
}

Gnuplotting & Gnuplotting::operator<<(const char * command){
	cmd(command);
	return *this;
}

void Gnuplotting::xrange(double xmin, double xmax){
	stringstream ss;
	ss << "set xrange [" << xmin << ":" << xmax << "]";
	string s = ss.str();
	cmd(s);
}

void Gnuplotting::yrange(double ymin, double ymax){
	stringstream ss;
	ss << "set yrange [" << ymin << ":" << ymax << "]";
	string s = ss.str();
	cmd(s);
}

void Gnuplotting::xlabel(string & x){
	string s = "set xlabel '" + x + "'";
	cmd(s);
}

void Gnuplotting::ylabel(string & y){
	string s = "set ylabel '" + y + "'";
	cmd(s);
}

void Gnuplotting::xlabel(const char * x){
	stringstream ss;
	ss << "set xlabel '" << x << "'";
	string s = ss.str();
	cmd(s);
}

void Gnuplotting::ylabel(const char * y){
	stringstream ss;
	ss << "set ylabel '" << y << "'";
	string s = ss.str();
	cmd(s);
}

void Gnuplotting::xystream(vector<double> & x, vector<double> & y){
	size_t N = x.size();
	cmd("plot '-' w lines");
	for(size_t i=0; i<N; i++){
		stringstream ss;
		ss << x[i] << " " << y[i];
		string str = ss.str();
		cmd(str);
	}
	cmd("e");
}

void Gnuplotting::xystream(size_t & N, arma::vec & x, arma::vec & y){
	cmd("plot '-' w lines");
	for(size_t i=0; i<N; i++){
		stringstream ss;
		ss << x(i) << " " << y(i);
		string str = ss.str();
		cmd(str);
	}
	cmd("e");
}

void Gnuplotting::xystream_replot(size_t & N, arma::vec & x, arma::vec & y, const char* title){
	stringstream s;
	s << "replot '-' title '" << title << "' w lines";
	string command_string = s.str();
	cout << command_string;
	cmd(command_string);
	for(size_t i=0; i<N; i++){
		stringstream ss;
		ss << x(i) << " " << y(i);
		string str = ss.str();
		cmd(str);
	}
	cmd("e");
}

void Gnuplotting::two_xystream(size_t &N1, vec &x1, vec &y1, const char* title1, size_t &N2, vec &x2, vec &y2, const char* title2){
	stringstream s;
	s << "plot '-' t '" << title1 << "' w lines, '-' t '" << title2 << "' w lines"; 
	string command1 = s.str();
	cmd(command1);
	for(size_t i=0; i < N1; i++){
		stringstream ss1;
		ss1 << x1(i) << " " << y1(i);
		string str1= ss1.str();
		cmd(str1);
	}
	cmd("e");
	for(size_t j=0; j<N2; j++){
		stringstream ss2;
		ss2 << x2(j) << " " << y2(j);
		string str2 = ss2.str();
		cmd(str2);
	}
	cmd("e");
}
