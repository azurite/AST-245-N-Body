#include <iostream>
#include <data.hpp>

Data::Data(const std::string &filename)
{
  // This is just the file format, note that we only need the first of the three numbers
  int numParticles, numGasParticles, numStarParticles;
  float curr;

  std::ifstream infile("data/" + filename);
  infile >> numParticles >> numGasParticles >> numStarParticles;

  size = numParticles;
  totalMass = 0.0;
  radius = 0.0;

  if(infile.is_open()) {
    // calculate the total mass of the system on the fly
    for(int i = 0; i < size; i++) {
      infile >> curr;
      mass.push_back(curr);
      totalMass += curr;
    }

    for(int i = 0; i < size; i++) {
      infile >> curr;
      x.push_back(curr);
    }

    for(int i = 0; i < size; i++) {
      infile >> curr;
      y.push_back(curr);
    }

    for(int i = 0; i < size; i++) {
      infile >> curr;
      z.push_back(curr);
    }

    // calculate the radius of our system. This is the largest radius of a particle from the center
    for(int i = 0; i < size; i++) {
      radius = std::max(radius, x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);
    }

    radius = std::sqrt(radius);

    for(int i = 0; i < size; i++) {
      infile >> curr;
      vx.push_back(curr);
    }

    for(int i = 0; i < size; i++) {
      infile >> curr;
      vy.push_back(curr);
    }

    for(int i = 0; i < size; i++) {
      infile >> curr;
      vz.push_back(curr);
    }

    for(int i = 0; i < size; i++) {
      infile >> curr;
      softening.push_back(curr);
    }

    for(int i = 0; i < size; i++) {
      infile >> curr;
      potential.push_back(curr);
    }
  }
  else {
    std::cout << "Could not open file" << std::endl;
  }

  std::cout << "Data successfully created" << std::endl;
}
