#include "functions.h"

pressure_function_t::pressure_function_t (
    const pressure_function_type_t type,
    const double C_rho): m_C_rho(C_rho), m_type(type)
{
  switch (type)
    {
    case pressure_function_type_t::linear:
      m_p = p_linear;
      m_d_p_d_rho = d_p_linear_d_rho;
      break;

    case pressure_function_type_t::exponential:
      m_p = p_gamma;
      m_d_p_d_rho = d_p_gamma_d_rho;
      break;

    default:
      break;
    }
}

double pressure_function_t::p (const double rho) const
{
  return m_p (rho, m_C_rho);
}

double pressure_function_t::d_p_d_rho (const double rho) const
{
  return m_d_p_d_rho (rho, m_C_rho);
}
