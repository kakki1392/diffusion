#include "gnuplotting.h"
#include <armadillo>

using namespace std;
using namespace arma;

int main(){


Gnuplotting plt;
double ymin = 0.0;
double ymax = 1.0;
plt.yrange(ymin,ymax);
plt.cmd("plot sin(x)");








	return 0;
}
