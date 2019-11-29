#include <iostream>
#include <data.hpp>
#include <particle.hpp>

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

std::vector<Particle> *Data::readFromFile(const std::string &filename)
{
  std::vector<float> mass;
  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;
  std::vector<float> vx;
  std::vector<float> vy;
  std::vector<float> vz;
  std::vector<float> softening;
  std::vector<float> potential;

  // This is just the file format, note that we only need the first of the three numbers
  int numParticles, numGasParticles, numStarParticles;
  float curr;

  std::ifstream infile("data/" + filename);
  infile >> numParticles >> numGasParticles >> numStarParticles;

  if(infile.is_open()) {
    // calculate the total mass of the system on the fly
    for(int i = 0; i < numParticles; i++) {
      infile >> curr;
      mass.push_back(curr);
    }

    for(int i = 0; i < numParticles; i++) {
      infile >> curr;
      x.push_back(curr);
    }

    for(int i = 0; i < numParticles; i++) {
      infile >> curr;
      y.push_back(curr);
    }

    for(int i = 0; i < numParticles; i++) {
      infile >> curr;
      z.push_back(curr);
    }

    for(int i = 0; i < numParticles; i++) {
      infile >> curr;
      vx.push_back(curr);
    }

    for(int i = 0; i < numParticles; i++) {
      infile >> curr;
      vy.push_back(curr);
    }

    for(int i = 0; i < numParticles; i++) {
      infile >> curr;
      vz.push_back(curr);
    }

    for(int i = 0; i < numParticles; i++) {
      infile >> curr;
      softening.push_back(curr);
    }

    for(int i = 0; i < numParticles; i++) {
      infile >> curr;
      potential.push_back(curr);
    }

    std::vector<Particle> *particles = new std::vector<Particle>();

    for(int i = 0; i < numParticles; i++) {
      Particle p(mass[i], x[i], y[i], z[i], vx[i], vy[i], vz[i], softening[i], potential[i]);
      particles->push_back(p);
    }

    return particles;
  }
  else {
    std::cout << "Could not open file" << std::endl;
    return NULL;
  }
}
