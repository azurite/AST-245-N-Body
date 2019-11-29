#include <iostream>
#include <Eigen/Dense>
#include <particle.hpp>

using Eigen::Vector3f;

Particle::Particle()
{
  mass = 0.0;
  pos << 0.0, 0.0, 0.0;
  velocity << 0.0, 0.0, 0.0;
  softening = 0.0;
  potential = 0.0;
}

Particle::Particle(float m, float x, float y, float z, float Vx, float Vy, float Vz, float s, float p)
{
  mass = m;
  pos << x, y, z;
  velocity << Vx, Vy, Vz;
  softening = s;
  potential = p;
}

float Particle::m()
{
  return mass;
}

void Particle::set_m(float m)
{
  mass = m;
}

float Particle::x()
{
  return pos(0);
}

void Particle::set_x(float x)
{
  pos(0) = x;
}

float Particle::y()
{
  return pos(1);
}

void Particle::set_y(float y)
{
  pos(1) = y;
}

float Particle::z()
{
  return pos(2);
}

void Particle::set_z(float z)
{
  pos(2) = z;
}

Vector3f Particle::r()
{
  return pos;
}

void Particle::set_r(Vector3f r)
{
  pos = r;
}

float Particle::Vx()
{
  return velocity(0);
}

void Particle::set_Vx(float Vx)
{
  velocity(0) = Vx;
}

float Particle::Vy()
{
  return velocity(1);
}

void Particle::set_Vy(float Vy)
{
  velocity(1) = Vy;
}

float Particle::Vz()
{
  return velocity(2);
}

void Particle::set_Vz(float Vz)
{
  velocity(2) = Vz;
}

Vector3f Particle::v()
{
  return velocity;
}

void Particle::set_v(Vector3f v)
{
  velocity = v;
}

float Particle::soft()
{
  return softening;
}

void Particle::set_soft(float s)
{
  softening = s;
}

float Particle::pot()
{
  return potential;
}

void Particle::set_pot(float p)
{
  potential = p;
}

std::string Particle::toString() const
{
  return "{ mass=" +
  std::to_string(mass) + ", x=" +
  std::to_string(pos(0)) + ", y=" +
  std::to_string(pos(1)) + ", z=" +
  std::to_string(pos(2)) + ", Vx=" +
  std::to_string(velocity(0)) + ", Vy=" +
  std::to_string(velocity(1)) + ", Vz=" +
  std::to_string(velocity(2)) + " }";
}

void Particle::print() const
{
  std::cout << toString() << std::endl;
}

std::ostream &operator<<(std::ostream &os, const Particle &p)
{
  os << p.toString() << std::endl;
  return os;
}
