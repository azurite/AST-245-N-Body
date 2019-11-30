#include <iostream>
#include <memory>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <mgl2/mgl.h>

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
  }

  return M / V;
}

float density_hernquist(float r)
{
  return totalMass / (2 * M_PI) * (scaleLength / r) * (1 / std::pow(r + scaleLength, 3));
}

void step1_calculate()
{
  // compute total mass and radius of the system
  for(Particle &p : particles) {
    totalMass += p.m();
    radius = std::max(p.radius2(), radius);
  }

  radius = std::sqrt(radius);
  scaleLength = radius * 0.45;

  std::cout << "First Task STEP 1 " << std::endl;
  std::cout << "------------------" << std::endl;
  std::cout << "  totalMass:      " << totalMass << std::endl;
  std::cout << "  radius:         " << radius << std::endl;
  std::cout << "  scaleLength:    " << scaleLength << std::endl;
  std::cout << "------------------" << std::endl;
}

void step1_plot()
{
  int numSteps = 1000;
  std::vector<float> hDensity;
  std::vector<float> nDensity;
  std::vector<float> rInput;

  for(float r = 0.001; r < radius; r += (radius / numSteps)) {
      float h_rho = density_hernquist(r);
      float n_rho = density_numeric(r);

      rInput.push_back(r);
      hDensity.push_back(h_rho);
      nDensity.push_back(n_rho);
  }

  mglData hData;
  hData.Set(hDensity.data(), hDensity.size());

  mglData nData;
  nData.Set(nDensity.data(), nDensity.size());

  mglData rData;
  rData.Set(rInput.data(), rInput.size());

  mglGraph gr;

  gr.SetRange('x', rData);
  gr.SetRange('y', hData);
  gr.SetFunc("", "lg(y)");
  gr.Adjust("y");
  gr.Axis();

  gr.Plot(hData, "b");
  gr.AddLegend("Hernquist", "b");

  gr.Plot(nData, "r");
  gr.AddLegend("Numeric", "r");

  gr.Legend();
  gr.WriteFrame("density_profiles.png");
}

void first_task_step1()
{
    step1_calculate();
    step1_plot();
}

void first_task_step2()
{

}

void first_task()
{
  first_task_step1();
  //first_task_step2();
}
