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


class pressure_function_t
{
  double (*m_p)(double, double) = nullptr;
  double (*m_d_p_d_rho) (double, double) = nullptr;
  double m_C_rho = 0.;

  pressure_function_type_t m_type = pressure_function_type_t::COUNT;

public:
  pressure_function_t (const pressure_function_type_t type, const double C_rho = 0.);
  pressure_function_t () = default;

  double p (const double rho) const;
  double d_p_d_rho (const double rho) const;

  pressure_function_type_t get_type () const { return m_type; }
  double get_C_rho () const { return m_C_rho; }
};

inline double rho_0_first (const double x)
{
  if (x < 4.5 || x > 5.5)
    return 1.;
  return 2.;
}
inline double u_0_first (const double /*x*/) { return 0.; }

inline double rho_0_second (const double /*x*/) { return 1.; }
inline double u_0_second (const double x)
{
  if (x < 4.5 || x > 5.5)
    return 0.;
  return 1.;
}

inline double zero_function (const double /*t*/, const double /*x*/) { return 0.; }

#endif // FUNCTIONS_H
