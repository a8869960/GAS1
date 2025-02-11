#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H

#include "functions.h"
#include "enums.h"

struct initial_conditions_t
{
  double (*m_f)(double, double, double, equation_of_state_t) = nullptr;
  double (*m_f0)(double, double) = nullptr;

  double m_mu = 0.;

  double (*m_rho)(double, double) = nullptr;
  double (*m_u)(double, double) = nullptr;

  initial_conditions_t () {};
  initial_conditions_t (const programm_mode_t type)
  {
    assert (type == programm_mode_t::debugging_test);

    m_f = f_debugging_test;
    m_f0 = f0_debugging_test;

    m_rho = rho_debugging_test;
    m_u = u_debugging_test;
  }

  double rho_0 (const double x) const { return m_rho (0., x); }
  double u_0 (const double x) const { return m_u (0., x); }
  double rho (const double t, const double x) const { return m_rho (t, x); }
  double u (const double t, const double x) const { return m_u (t, x); }
  double f (const double t, const double x) const  { return m_f (t, x, m_mu, m_p); }
  double f0 (const double t, const double x) const { return m_f0 (t, x); }

  equation_of_state_t m_p;
};

#endif // INITIAL_CONDITIONS_H
