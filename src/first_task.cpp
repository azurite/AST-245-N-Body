#include <iostream>
#include <memory>
#include <algorithm>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <mgl2/mgl.h>
#include <Eigen/Dense>

#include <data.hpp>
#include <particle.hpp>
#include <gravitysolvers.hpp>
#include <first_task.hpp>

std::unique_ptr<std::vector<Particle>> p(Data::readFromFile("data.ascii"));
std::vector<Particle> particles = *p;

float pMass = 0.0;
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

float density_hernquist(float r)
{
  return (totalMass / (2 * M_PI)) * (scaleLength / r) * (1 / std::pow(r + scaleLength, 3));
}

float force_hernquist(float r)
{
  // Assumption: G = 1
  return -totalMass * pMass / ((r + scaleLength) * (r + scaleLength));
}

void calculate_constants()
{
  pMass = particles[0].m(); // all particles have the same mass

  for(Particle &p : particles) {
    totalMass += p.m();
    radius = std::max(p.radius2(), radius);
    r0 = std::min(p.radius2(), r0);
  }

  r0 = std::sqrt(r0) + std::numeric_limits<float>::epsilon();
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
  std::cout << "  pMass:          " << pMass << std::endl;
  std::cout << "  totalMass:      " << totalMass << std::endl;
  std::cout << "  radius:         " << radius << std::endl;
  std::cout << "  r0:             " << r0 << std::endl;
  std::cout << "  rhm:            " << rhm << std::endl;
  std::cout << "  scaleLength:    " << scaleLength << std::endl;
  std::cout << "  softening:      " << epsilon << std::endl;
  std::cout << "  t_relax:        " << t_relax << std::endl;
  std::cout << "------------------" << std::endl;
}

float compute_relaxation()
{
  // Assumption: G = 1
  float N = particles.size();
  float vc2 = rhm * std::abs(force_hernquist(rhm));
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
      float r0 = r;
      float r1 = r + dr;
      float MShell = Mass(r1, true) - Mass(r0, true);
      float VShell = 4.0/3.0*M_PI*(r1*r1*r1 - r0*r0*r0);

      // we expand the radius of the shell until we find nonzero mass inside of it
      // in order to yield visually more meaningful results on a log plot
      while(MShell <= std::numeric_limits<float>::epsilon()) {
          r1 += dr;
          MShell = Mass(r1, true) - Mass(r0, true);
      }

      float n_rho = MShell / VShell;
      float numParticles = MShell / pMass;

      // p = n*m/v, err = sqrt(n) =>
      // p_err = sqrt(n)/n*p = sqrt(n)*n*m/(n*v) = sqrt(n)*m/v and adjust error to log10 scale
      float rho_error = std::sqrt(numParticles) * pMass / VShell / std::log(10);

      hDensity.push_back(density_hernquist((r0 + r1) / 2));
      nDensity.push_back(n_rho);
      errors.push_back(rho_error);

      rInput.push_back(r);
      r = r1 - dr;
  }

  mglData hData;
  hData.Set(hDensity.data(), hDensity.size());

  mglData nData;
  nData.Set(nDensity.data(), nDensity.size());

  mglData rData;
  rData.Set(rInput.data(), rInput.size());

  mglData eData;
  eData.Set(errors.data(), errors.size());

  mglGraph gr(0, 1200, 800);

  float outMin = std::min(hData.Minimal(), nData.Minimal());
  float outMax = std::max(hData.Maximal(), nData.Maximal());

  gr.SetRange('x', rData);
  gr.SetRange('y', outMin, outMax);
  gr.SetCoor(mglLogY);
  gr.Axis();

  gr.Label('x', "Radius [l]", 0);
  gr.Label('y', "Density [m]/[l]^3", 0);

  gr.Plot(rData, hData, "b");
  gr.AddLegend("Hernquist", "b");

  gr.Plot(rData, nData, "r .");
  gr.AddLegend("Numeric", "r .");

  gr.Error(rData, nData, eData, "qo");
  gr.AddLegend("Poissonian Error", "qo");

  gr.Legend();
  gr.WritePNG("density_profiles.png");
}

void step2()
{
  r0 = 0.005;
  int numSteps = 100;

  std::vector<float> softenings;
  std::vector<float> dAnalytic;
  std::vector<float> rInput;

  // creates evenly spaced intervals on a log scale on the interval [r0, radius]
  // drLinToLog(i) gives the start of the ith interval on [r0, radius]
  auto drLinToLog = [&](int i) {
    return r0 * std::pow(radius / r0, (float)i / (float)numSteps);
  };

  // transforms the numerical data into a dimension desirable for plotting
  // here we want to plot the magnitude of the force and add a small nonzero
  // number to the force to avoid problems with 0 on log scales since log(0) is undefined
  auto plotFit = [](float x) {
    return std::abs(x) + std::numeric_limits<float>::epsilon();
  };

  // calculate the analytical force in the hernquist model
  for(int i = 0; i <= numSteps; i++) {
    float r = drLinToLog(i);
    dAnalytic.push_back(plotFit(force_hernquist(r)));
    rInput.push_back(r);
  }

  std::unique_ptr<Gravitysolver::Direct> solver(new Gravitysolver::Direct());
  std::vector<mglData> plotData;

  for(int i = 1; i <= 7; i++) {
    solver->readData("data/direct-nbody-" + std::to_string(i) + ".txt");
    softenings.push_back(solver->softening());

    std::vector<float> dNumeric;
    Eigen::VectorXf fn = solver->data().row(10);

    // calculate the average gravitational force in a shell
    for(int i = 0; i <= numSteps; i++) {
      float r = drLinToLog(i);
      float r1 = drLinToLog(i + 1);
      float dr = r1 - r;

      float f = 0.0;
      int numParticles = 0;

      for(int i = 0; i < particles.size(); i++) {
        Particle p = particles[i];

        if(r <= p.radius() && p.radius() < r1) {
          f += fn(i);
          numParticles++;
        }
      }

      f = (numParticles == 0 ? .0 : f / numParticles);
      dNumeric.push_back(plotFit(f));
    }

    mglData cData;
    cData.Set(dNumeric.data(), dNumeric.size());
    plotData.push_back(cData);
  }

  mglData aData;
  aData.Set(dAnalytic.data(), dAnalytic.size());

  mglData rData;
  rData.Set(rInput.data(), rInput.size());

  mglGraph gr(0, 1200, 800);

  float outMin;
  float outMax;

  for(int i = 0; i < plotData.size(); i++) {
    outMin = std::max(std::min(plotData[i].Minimal(), aData.Minimal()), 50.0);
    outMax = std::max(plotData[i].Maximal(), aData.Maximal());
  }

  gr.SetRange('x', rData);
  gr.SetRange('y', outMin, outMax);

  gr.SetFontSize(2);
  gr.SetCoor(mglLogLog);
  gr.Axis();

  gr.Label('x', "Radius [l]", 0);
  gr.Label('y', "Force [m]^2[l]^{-2}", 0);

  gr.Plot(rData, aData, "b");
  gr.AddLegend("Analytic", "b");

  // colors for plotting
  const char *opt[7] = {"r +", "c +", "m +", "h +", "l +", "n +", "q +"};

  for(int i = 0; i < plotData.size(); i++) {
    std::stringstream ss;
    ss << "\\epsilon = " << std::setprecision(4) << softenings[i] << " [l]";

    gr.Plot(rData, plotData[i], opt[i]);
    gr.AddLegend(ss.str().c_str(), opt[i]);
  }

  gr.Legend();
  gr.WritePNG("forces.png");
}

void first_task()
{
  calculate_constants();
  step1();
  step2();
}
