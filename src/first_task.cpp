#include <iostream>
#include <memory>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <data.hpp>
#include <particle.hpp>
#include <first_task.hpp>

std::unique_ptr<std::vector<Particle>> p(Data::readFromFile("data.ascii"));
std::vector<Particle> particles = *p;

float totalMass = 0.0;
float radius = 0.0;
float scaleLength = 0.0;

float density_numeric(float r)
{
  float M = 0.0;
  float V = 4/3*M_PI*r*r*r;

  for(int i = 0; i < particles.size(); i++) {
    Particle p = particles[i];

    if(p.radius2() <= r * r) {
      M += p.m();
    }
    else {
      break;
    }
  }

  return M / V;
}

float density_hernquist(float r)
{
  return totalMass / (2 * M_PI) * (scaleLength / r) * (1 / std::pow(r + scaleLength, 3));
}

// https://astronomy.swin.edu.au/cms/astro/cosmos/S/Scale+Length
// TODO: This is probably not correct as the scale length is neglegile compared to the radius
// try the 0.45 * radius from the paper
float calc_scale_length()
{
  float r = 0.05;
  float rho_center = density_numeric(r);
  float rho_curr = rho_center;

  while(rho_center < M_E * rho_curr) {
    r += 0.05;
    rho_curr = density_numeric(r);
  }

  return r;
}

void first_task()
{
  // sort particles from the inner most to the outer most
  std::sort(particles.begin(), particles.end(), [](Particle &a, Particle &b){
    //return a.x()*a.x() + a.y()*a.y() + a.z()*a.z() < b.x()*b.x() + b.y()*b.y() + b.z()*b.z();
    return a.radius2() < b.radius2();
  });

  // compute total mass and radius of the system
  for(Particle &p : particles) {
    totalMass += p.m();
  }

  radius = particles.back().radius();
  scaleLength = calc_scale_length();

  std::cout << "totalMass: " << totalMass << std::endl;
  std::cout << "radius: " << radius << std::endl;
  std::cout << "scaleLength: " << scaleLength << std::endl;
}
