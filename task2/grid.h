#ifndef GRID_H
#define GRID_H

#define EPS 1e-16

#include <cassert>
#include "math.h"

class grid_t
{
private:
  double m_tau;
  double m_h;

  double m_t = 0.;
  double m_T = 1.;
  double m_x = 0.;
  double m_X = 10.;

public:
  grid_t () {}
  grid_t (const double tau, const double h): m_tau(tau), m_h(h) {}
  grid_t (const double tau, const double h, const double T, const double X): m_tau(tau), m_h(h), m_T(T), m_X(X) {}

  double get_tau () const { return m_tau; }
  double get_h () const { return m_h; }

  double get_x () const { return m_x; }
  double get_X () const { return m_X; }
  double get_t () const { return m_t; }
  double get_T () const { return m_T; }

  unsigned int get_M () const { return (unsigned int)(m_X - m_x) / m_h + 1; }
  unsigned int get_N () const { return (unsigned int)(m_T - m_t) / m_tau; }

  double t (const unsigned int n) const { return m_t + n * m_tau; }
  double m (const unsigned int n) const { return m_x + n * m_h; }
};


#endif // GRID_H
