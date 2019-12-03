#include <iostream>
#include <memory>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <mgl2/mgl.h>
#include <Eigen/Dense>

#include <data.hpp>
#include <particle.hpp>
#include <first_task.hpp>

using Eigen::Vector3f;

void writeNBody(std::vector<float> *a, float eps, const std::string &filename)
{
  std::ofstream outfile("data/" + filename);

  if(outfile.is_open()) {
    outfile << a->size() << " " << eps << "\n";

    for(int i = 0; i < a->size(); i++) {
      outfile << a->at(i) << "\n";
    }

    outfile.close();
  }
  else {
    std::cout << "could not write to file" << std::endl;
  }
}

std::vector<float> *readNBody(const std::string &filename)
{
  std::ifstream infile("data/" + filename);

  if(infile.is_open()) {
    int size;
    infile >> size;

    std::vector<float> *a = new std::vector<float>(size + 1, 0.0);
    infile >> (*a)[0];

    for(int i = 1; i < size + 1; i++) {
      infile >> (*a)[i];
    }

    return a;
  }
  else {
    std::cout << "could not read from file" << std::endl;
    return NULL;
  }
}

std::unique_ptr<std::vector<Particle>> p(Data::readFromFile("data.ascii"));
std::vector<Particle> particles = *p;

float totalMass = 0.0;
float radius = 0.0;
float scaleLength = 0.0;
float epsilon = 0.0; // softening
float r0 = 0.1; // center of system to avoid problems with dividing by 0

// because computing the N-Boby forces is an n^2 algorithm it can take quite long to run depending
// on the hardware. The blockSize variable specifies what fraction of particles
// is going to be used. F.ex a blockSize of 10 means we only take every 10th element
// if blockSize = 1 this corresponds to taking all particles into consideration
int blockSize = 1;

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

void calculate_constants()
{
  for(Particle &p : particles) {
    totalMass += p.m();
    radius = std::max(p.radius2(), radius);
  }

  radius = std::sqrt(radius);
  scaleLength = radius * 0.45;

  // https://en.wikipedia.org/wiki/Mean_inter-particle_distance
  // after plugging in n = totalMass / (4/3*PI*r^3) in the formula most terms cancel out
  epsilon = radius / std::pow(totalMass, 1.0/3.0);

  std::cout << "First Task        " << std::endl;
  std::cout << "------------------" << std::endl;
  std::cout << "  totalMass:      " << totalMass << std::endl;
  std::cout << "  radius:         " << radius << std::endl;
  std::cout << "  scaleLength:    " << scaleLength << std::endl;
  std::cout << "  softening:      " << epsilon << std::endl;
  std::cout << "------------------" << std::endl;
}

std::vector<float> *force_n_body(const std::string &fname)
{
  std::vector<float> *a = readNBody(fname);

  if(a != NULL) {
    epsilon = a->front();
    a->erase(a->begin());

    std::cout << "using n-body forces from file" << std::endl;
    return a;
  }

  a = new std::vector<float>(particles.size(), 0.0);

  for(int i = 0; i <= particles.size() - blockSize; i += blockSize) {
    Vector3f ai(.0, .0, .0);

    for(int j = 0; j < particles.size(); j++) {
      if(i != j) {
        Vector3f rij = particles[j].r() - particles[i].r();
        float x = (rij.squaredNorm() + epsilon * epsilon);

        // Assumption: G = 1
        ai += particles[j].m() * rij / (x * std::sqrt(x));
      }
    }

    // project the force vector onto the normalized sphere normal
    (*a)[i] = ai.dot(particles[i].r() / particles[i].radius());
  }

  writeNBody(a, epsilon, "nbody_a.ascii");

  return a;
}

float force_analytic(float r)
{
  float M = 0.0;
  float r2 = r * r;

  for(int i = 0; i < particles.size(); i++) {
    if(particles[i].radius2() <= r2) {
      M += particles[i].m();
    }
  }

  // Assumption: G = 1
  float x = (r2 + epsilon * epsilon);
  return (-M * r) / (x * std::sqrt(x));
}

void step1()
{
  int numSteps = 1000;
  std::vector<float> hDensity;
  std::vector<float> nDensity;
  std::vector<float> rInput;

  for(float r = r0; r <= radius; r += (radius / numSteps)) {
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

  mglGraph gr(0, 1200, 800);

  float outMin = std::min(hData.Minimal(), nData.Minimal());
  float outMax = std::max(hData.Maximal(), nData.Maximal());

  gr.SetRange('x', rData);
  gr.SetRange('y', outMin, outMax);
  gr.SetFunc("", "lg(y)");
  gr.Adjust("y");
  gr.Axis();

  gr.Plot(hData, "b");
  gr.AddLegend("Hernquist", "b");

  gr.Plot(nData, "r");
  gr.AddLegend("Numeric", "r");

  gr.Legend();
  gr.WritePNG("density_profiles.png");
}

void step2()
{
  std::unique_ptr<std::vector<float>> ptr(force_n_body("nbody_a.ascii"));
  std::vector<float> a_numeric = *ptr;

  int numSteps = 1000;
  float dr = radius / numSteps;
  std::vector<float> dAnalytic;
  std::vector<float> dNumeric;
  std::vector<float> rInput;

  for(float r = r0; r <= radius; r += dr) {
    float a_curr = 0.0;
    int numParticles = 0;

    for(int i = 0; i < particles.size(); i++) {
      Particle p = particles[i];

      if((r - dr/2) < p.radius() && p.radius() <= (r + dr/2)) {
        a_curr += a_numeric[i];
        numParticles++;
      }
    }

    // take the average force felt by all particles in the shell
    a_curr = numParticles == 0 ? 0.0 : a_curr / numParticles;

    dAnalytic.push_back(force_analytic(r));
    dNumeric.push_back(a_curr);
    rInput.push_back(r);
  }

  mglData aData;
  aData.Set(dAnalytic.data(), dAnalytic.size());

  mglData nData;
  nData.Set(dNumeric.data(), dNumeric.size());

  mglData rData;
  rData.Set(rInput.data(), rInput.size());

  mglGraph gr(0, 1200, 800);

  float outMin = std::min(aData.Minimal(), nData.Minimal());
  float outMax = std::max(aData.Maximal(), nData.Maximal());

  gr.SetRange('x', rData);
  gr.SetRange('y', outMin, outMax);
  gr.SetFunc("lg(x)", "");
  gr.Adjust("y");
  gr.Axis();

  gr.Plot(aData, "b");
  gr.AddLegend("Analytic", "b");

  gr.Plot(nData, "r");
  gr.AddLegend("Numeric", "r");

  gr.Legend();
  gr.WritePNG("forces.png");
}

void first_task()
{
  calculate_constants();
  step1();
  step2();
}
