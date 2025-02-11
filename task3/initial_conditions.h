#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H

#include "functions.h"
#include "enums.h"

struct initial_conditions_t
{
  double (*m_f)(double, double) = nullptr;
  double (*m_f0)(double, double) = nullptr;

  double m_mu = 0.;
  unsigned int k;

  double (*m_rho)(double, double) = nullptr;
  double (*m_u)(double, double) = nullptr;

  double (*m_rho_0)(double, unsigned int) = nullptr;         // rho_0 (x)
  double (*m_u_0)(double, unsigned int) = nullptr;           // u_0 (x)

  initial_conditions_t () {};
  initial_conditions_t (const programm_mode_t mode)
  {
    m_f = zero_function;
    m_f0 = zero_function;

    if (mode == programm_mode_t::second)
      {
        m_rho_0 = rho_0_second;
        m_u_0 = u_0_second;
      }
    else
      {
        m_rho_0 = rho_0_first;
        m_u_0 = u_0_first;
      }
  }

  double rho_0 (const double x) const { return m_rho_0 (x, k); }
  double u_0 (const double x) const { return m_u_0 (x, k); }
  double rho (const double t, const double x) const { return m_rho (t, x); }
  double u (const double t, const double x) const { return m_u (t, x); }
  double f (const double t, const double x) const  { return m_f (t, x); }
  double f0 (const double t, const double x) const { return m_f0 (t, x); }

  pressure_function_t m_p;
};

#endif // INITIAL_CONDITIONS_H
