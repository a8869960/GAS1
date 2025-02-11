#include "functions.h"

equation_of_state_t::equation_of_state_t (
    const equation_of_state_type_t type,
    const double C_rho): m_C_rho(C_rho), m_type(type)
{
  switch (type)
    {
    case equation_of_state_type_t::linear:
      m_p = p_linear;
      m_d_p_d_rho = d_p_linear_d_rho;
      break;

    case equation_of_state_type_t::exponential:
      m_p = p_gamma;
      m_d_p_d_rho = d_p_gamma_d_rho;
      break;

    default:
      break;
    }
}

double equation_of_state_t::p (const double rho) const
{
  return m_p (rho, m_C_rho);
}

double equation_of_state_t::d_p_d_rho (const double rho) const
{
  return m_d_p_d_rho (rho, m_C_rho);
}
