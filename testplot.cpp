#include "gnuplotting.h"
#include <unistd.h>
#include <cstdlib>
#include <armadillo>

using namespace std;
using namespace arma;

int main(){


Gnuplotting plt;
double ymin = -1.1;
double ymax = 1.1;
plt.yrange(ymin,ymax);
plt.xrange(-0.1,1.1);

size_t N = 10;
vec v1 = randu<vec>(N);
vec x1 = randu<vec>(N);

vec v2 = randu<vec>(N);
vec x2 = randu<vec>(N);

plt.two_xystream(N,x1,v1,"t1",N,x2,v2,"t2");
v1 = 2.0*v1;
sleep(2);
plt.two_xystream(N,x1,v1,"t1",N,x2,v2,"t2");
	return 0;
}
