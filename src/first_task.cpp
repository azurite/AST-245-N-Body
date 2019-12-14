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
float t_relax = 0.0;
float r0 = std::numeric_limits<float>::max(); // center of system to avoid problems with dividing by 0
float rhm = 0.0; // half mass radius

// because computing the N-Boby forces is an n^2 algorithm it can take quite long to run depending
// on the hardware. The blockSize variable specifies what fraction of particles
// is going to be used. F.ex a blockSize of 10 means we only take every 10th element
// if blockSize = 1 this corresponds to taking all particles into consideration
int blockSize = 1;

float Mass(float r, bool strictlyLess = false)
{
  float M = 0.0;

  for(Particle &p : particles) {
    float r2 = p.radius2();

    if((!strictlyLess && r2 <= r*r) || (strictlyLess && r2 < r*r)) {
      M += p.m();
    }
  }

  return M;
}

float density_in_shell(float r, float dr)
{
  float m = .0;
  float r1 = r + dr;

  for(Particle &p : particles) {
    if(r <= p.radius() && p.radius() < r + dr) {
      m += p.m();
    }
  }

  return m / (4.0/3.0*M_PI*(r1*r1*r1 - r*r*r));
}

float density_hernquist(float r)
{
  return totalMass / (2 * M_PI) * (scaleLength / r) * (1 / std::pow(r + scaleLength, 3));
}

float force_analytic(float r)
{
  // Assumption: G = 1
  return -Mass(r) / (r * r);
}

void calculate_constants()
{
  for(Particle &p : particles) {
    totalMass += p.m();
    radius = std::max(p.radius2(), radius);
    r0 = std::min(p.radius2(), r0);
  }

  r0 = std::sqrt(r0);
  radius = std::sqrt(radius);

  float dr = radius / 100000;
  for(float r = r0; r <= radius; r += dr) {
    if(Mass(r) >= totalMass * 0.5) {
      rhm = r;
      break;
    }
  }

  scaleLength = rhm / (1 + std::sqrt(2));

  // https://en.wikipedia.org/wiki/Mean_inter-particle_distance
  // after plugging in n = totalMass / (4/3*PI*r^3) in the formula most terms cancel out
  epsilon = radius / std::pow(totalMass, 1.0/3.0);
  t_relax = compute_relaxation();

  std::cout << "First Task        " << std::endl;
  std::cout << "------------------" << std::endl;
  std::cout << "  totalMass:      " << totalMass << std::endl;
  std::cout << "  radius:         " << radius << std::endl;
  std::cout << "  r0:             " << r0 << std::endl;
  std::cout << "  rhm:            " << rhm << std::endl;
  std::cout << "  scaleLength:    " << scaleLength << std::endl;
  std::cout << "  softening:      " << epsilon << std::endl;
  std::cout << "  t_relax:        " << t_relax << std::endl;
  std::cout << "------------------" << std::endl;
}

std::vector<float> *force_n_body(const std::string &fname)
{
  std::vector<float> *a = readNBody(fname);
  bool usingFile = false;

  if(a != NULL) {
    usingFile = true;
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

  if(!usingFile) {
    writeNBody(a, epsilon, "nbody_a.ascii");
  }

  return a;
}

float compute_relaxation()
{
  // Assumption: G = 1
  float N = particles.size();
  float vc2 = rhm * std::abs(force_analytic(rhm));
  float t_cross = radius / std::sqrt(vc2);
  float t_relax = N / (8 * std::log(N)) * t_cross;

  return t_relax;
}

void step1()
{
  int numSteps = 100;
  float dr = radius / numSteps;
  std::vector<float> hDensity;
  std::vector<float> nDensity;
  std::vector<float> rInput;
  std::vector<float> errors;

  for(float r = r0; r < radius; r += dr) {
      float h_rho = density_hernquist(r + dr/2);
      float n_rho = density_in_shell(r, dr) + std::numeric_limits<float>::epsilon();

      float r1 = r + dr;
      float expectedDensity = (Mass(r1, true) - Mass(r, true)) / (4.0/3.0*M_PI*(r1*r1*r1 - r*r*r));

      rInput.push_back(r);
      hDensity.push_back(h_rho);
      nDensity.push_back(n_rho);
      errors.push_back(std::sqrt(expectedDensity) + std::numeric_limits<float>::epsilon());
  }

  mglData hData;
  hData.Set(hDensity.data(), hDensity.size());

  mglData nData;
  nData.Set(nDensity.data(), nDensity.size());

  mglData rData;
  rData.Set(rInput.data(), rInput.size());

  mglData eData;
  eData.Set(errors.data(), errors.size());

  mglGraph grDensity(0, 1200, 800);

  float outMin = std::min(hData.Minimal(), nData.Minimal());
  float outMax = std::max(hData.Maximal(), nData.Maximal());

  grDensity.SetRange('x', rData);
  grDensity.SetRange('y', outMin, outMax);
  grDensity.SetCoor(mglLogY);
  grDensity.Axis();

  grDensity.Plot(hData, "b");
  grDensity.AddLegend("Hernquist", "b");

  grDensity.Plot(nData, "r +");
  grDensity.AddLegend("Numeric", "r +");

  //grDensity.Error(rData, nData, eData);

  grDensity.Legend();
  grDensity.WritePNG("density_profiles.png");

  mglGraph grError(0, 1200, 800);

  grError.SetRange('x', rData);
  grError.SetRange('y', eData);
  grError.SetCoor(mglLogY);
  grError.Axis();

  grError.Plot(eData, "r +");
  grError.AddLegend("Poissonian Error", "r +");

  grError.Legend();
  grError.WritePNG("poissonain_error.png");
}

void step2()
{
  std::unique_ptr<std::vector<float>> ptr(force_n_body("nbody_a.ascii"));
  std::vector<float> a_numeric = *ptr;

  int numSteps = 2000;
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

  // compute the relaxation time scale for different softenings
  float epsilon0 = epsilon;
  std::vector<float> dRelax;
  std::vector<float> dSoft;

  for(float eps = epsilon0 / 8; eps <= epsilon0 * 8; eps *= 2) {
    epsilon = eps;
    dRelax.push_back(compute_relaxation());
    dSoft.push_back(epsilon);
  }

  epsilon = epsilon0;

  mglData relaxData;
  relaxData.Set(dRelax.data(), dRelax.size());

  mglData softData;
  softData.Set(dSoft.data(), dSoft.size());

  mglGraph gr_relax(0, 1200, 800);

  gr_relax.SetRange('x', softData);
  gr_relax.SetRange('y', relaxData);
  gr_relax.Axis();

  gr_relax.Plot(relaxData, "b");
  gr_relax.AddLegend("t_{relax}", "b");

  gr_relax.Legend();
  gr_relax.WritePNG("relaxation.png");
}

void first_task()
{
  calculate_constants();
  step1();
  //step2();
}
