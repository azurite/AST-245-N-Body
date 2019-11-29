#include <iostream>
#include <Eigen/Dense>

using Eigen::Vector3f;

#ifndef PARTICLE_H
#define PARTICLE_H

class Particle
{
  float mass; // Mass
  Vector3f pos; // Position Vector (x, y, z)
  Vector3f velocity; // Velocity Vector (vx, vy, vz)
  float softening;
  float potential;

public:
  Particle();
  Particle(float m, float x, float y, float z, float Vx, float Vy, float Vz, float softening, float potential);

  float m();
  void set_m(float m);
  float x();
  void set_x(float x);
  float y();
  void set_y(float y);
  float z();
  void set_z(float z);
  Vector3f r();
  void set_r(Vector3f r);

  float Vx();
  void set_Vx(float Vx);
  float Vy();
  void set_Vy(float Vy);
  float Vz();
  void set_Vz(float Vz);
  Vector3f v();
  void set_v(Vector3f v);

  float soft();
  void set_soft(float s);
  float pot();
  void set_pot(float p);

  std::string toString() const;
  void print() const;

  friend std::ostream &operator<<(std::ostream &os, const Particle &p);
};

#endif // PARTICLE_H
