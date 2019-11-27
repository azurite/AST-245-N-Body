#include <fstream>
#include <vector>
#include "particle.hpp"

#ifndef DATA_H
#define DATA_H

/*
* The data file is made of a header and sequential arrays(each array as a single column data entry), as follows: Header
* Number of particles(N), Number of gas particles(N*=0 there is actually no gas particle), Number of star particles
* Arrays(running in i=1,...,N)
* Masses[i]
* x[i]
* y[i]
* z[i]
* Vx[i]
* Vy[i]
* Vz[i]
* softening[i]
* potential[i]
*/

class Data
{
  int size;
  std::vector<float> mass;
  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;
  std::vector<float> vx;
  std::vector<float> vy;
  std::vector<float> vz;
  std::vector<float> softening;
  std::vector<float> potential;

  float totalMass;
  float radius;

public:
  Data(const std::string &filename);
};

#endif // DATA_H
