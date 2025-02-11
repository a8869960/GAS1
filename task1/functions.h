#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "math.h"
#include "enums.h"

#include <cassert>

#define gamma 1.4

inline double p_linear (const double rho, const double C_rho) { return C_rho * rho; }
inline double d_p_linear_d_rho (const double /*rho*/, const double C_rho) { return C_rho; }
inline double p_gamma (const double rho, const double) { if (rho < 0) return 1.; return pow (rho, gamma); }
inline double d_p_gamma_d_rho (const double rho, const double) { return pow (rho, gamma - 1.); }


class equation_of_state_t
{
  double (*m_p)(double, double) = nullptr;
  double (*m_d_p_d_rho) (double, double) = nullptr;
  double m_C_rho = 0.;

  equation_of_state_type_t m_type = equation_of_state_type_t::COUNT;

public:
  equation_of_state_t (const equation_of_state_type_t type, const double C_rho = 0.);
  equation_of_state_t () = default;

  double p (const double rho) const;
  double d_p_d_rho (const double rho) const;

  equation_of_state_type_t get_type () const { return m_type; }
  double get_C_rho () const { return m_C_rho; }
};

// Debugging test
inline double rho_debugging_test (const double t, const double x)
{
  return exp (t) * (cos (M_PI * x) + 1.5);
}
inline double u_debugging_test (const double t, const double x)
{
  return cos (M_PI * t) * sin (M_PI * x);
}

inline double f0_debugging_test (const double t, const double x)
{
  return rho_debugging_test (t, x) + rho_debugging_test (t, x) * M_PI * cos (M_PI * t) * cos (M_PI * x)
      - exp (t) * M_PI * sin (M_PI * x) * u_debugging_test (t, x);
}

inline double f_debugging_test_linear (
    const double t,
    const double x,
    const double mu,
    const equation_of_state_t p)
{
  assert (p.get_type () == equation_of_state_type_t::linear);
  const double numerator = -p.get_C_rho () * exp (t) * M_PI * sin (M_PI * x)
      + mu * M_PI * M_PI * u_debugging_test (t, x)
      + M_PI * cos (M_PI * t) * cos (M_PI * x) * rho_debugging_test (t, x) * u_debugging_test (t, x)
      - exp (t) * M_PI * (1.5 + cos (M_PI * x)) * sin (M_PI * t) * sin (M_PI * x);
  return numerator / rho_debugging_test (t, x);
}

inline double f_debugging_test_expotential (
    const double t,
    const double x,
    const double mu,
    const equation_of_state_t p)
{
  assert (p.get_type () == equation_of_state_type_t::exponential);
  const double numerator = mu * M_PI * M_PI * u_debugging_test (t, x)
      + M_PI * u_debugging_test (t, x) * rho_debugging_test (t, x) * cos (M_PI * t) * cos (M_PI * x)
      - exp (t) * M_PI * gamma * pow (rho_debugging_test (t, x), gamma - 1) * sin (M_PI * x)
      - M_PI * rho_debugging_test (t, x) * sin (M_PI * t) * sin (M_PI * x);
  return numerator / rho_debugging_test (t, x);
}

inline double f_debugging_test (
    const double t,
    const double x,
    const double mu,
    const equation_of_state_t p)
{
  if (p.get_type () == equation_of_state_type_t::linear)
    return f_debugging_test_linear (t, x, mu, p);
  else
    return f_debugging_test_expotential (t, x, mu, p);
}

/*
inline double f0_debugging_test (const double t, const double x)
{
  const double d_rho_d_t = exp (t) * (cos (3. * M_PI * x) + 1.5);
  const double d_rho_u_d_x = exp (t) * cos (2. * M_PI * t)
      * (3. * M_PI * cos (7. * M_PI * x) + M_PI * cos (3. * M_PI * x) * cos (4. * M_PI * x));
  return d_rho_d_t + d_rho_u_d_x;
}
inline double f_debugging_test (
    const double t,
    const double x,
    const double mu,
    const equation_of_state_t p)
{
  const double rho = rho_debugging_test (t, x);
  const double u = u_debugging_test (t, x);
  const double d_u_d_t = - 2. * M_PI * sin (2. * M_PI * t) * sin (4. * M_PI * x);
  const double d_u_d_x = 4. * M_PI * cos (2. * M_PI * t) * cos (4. * M_PI * x);
  const double d_rho_d_x = - 3. * M_PI * exp (t) * sin (3. * M_PI * x);
  const double d2_u_d_x2 = - 16. * M_PI * M_PI * cos (2. * M_PI * t) * sin (4. * M_PI * x);
  return rho * d_u_d_t + rho * u * d_u_d_x + d_rho_d_x * p.d_p_d_rho (rho) - mu * d2_u_d_x2;
}
*/

#endif // FUNCTIONS_H
