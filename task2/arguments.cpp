#include "arguments.h"

error_io arguments_t::is_mu_correct () const
{
  if (mu < 0. + EPS || mu > 1 + EPS)
    return error_io::incorrect_arguments;

  // Supported only mu = 1e-1, 1e-2, 1e-3
  if (!(fabs (mu - 1e-1) < EPS || fabs (mu - 1e-2) < EPS || fabs (mu - 1e-3) < EPS))
    return error_io::unsupported_mu;

  return error_io::success;
}

error_io arguments_t::is_C_rho_correct () const
{
  if (C_rho < 0. - EPS || C_rho > 1e3 + EPS)
    return error_io::incorrect_arguments;

  // Supported only C_rho = 1, 10, 100
  /// If C_rho == 0 use expotential equation of state
  if (!(fabs (C_rho) < EPS || fabs (C_rho - 1.) < EPS || fabs (C_rho - 10.) < EPS || fabs (C_rho - 100.) < EPS))
    return error_io::unsupported_C_rho;

  return error_io::success;
}

bool arguments_t::is_correct () const
{
  error_io mu_status = is_mu_correct ();
  error_io C_rho_status = is_C_rho_correct ();

  print_error_io (mu_status);
  print_error_io (C_rho_status);

  const bool mu_correct = (mu_status == error_io::success || mu_status == error_io::unsupported_mu);
  const bool C_rho_correct = (C_rho_status == error_io::success || C_rho_status == error_io::unsupported_C_rho);

  return (mu_correct && C_rho_correct && tau > 0. && h > 0.);
}

error_io arguments_t::read (const int argc, char *argv[])
{
  if (argc != 10)
    return error_io::incorrect_arguments;

  int mode_task;
  if (   sscanf (argv[1], "%lf", &mu) != 1
      || sscanf (argv[2], "%lf", &C_rho) != 1
      || sscanf (argv[3], "%lf", &tau) != 1
      || sscanf (argv[4], "%lf", &h) != 1
      || sscanf (argv[5], "%d", &mode_task) != 1
      || sscanf (argv[8], "%u", &n_st) != 1
      || sscanf (argv[9], "%lf", &eps) != 1)
    return error_io::error_read;

  programm_mode_t programm_task = (mode_task == 1) ? programm_mode_t::first : programm_mode_t::second;
  const std::string p_name = argv[6];
  if (p_name == "linear" || p_name == "l")
    p = pressure_function_t (pressure_function_type_t::linear, C_rho);
  else if (p_name == "expotential" || p_name == "exp" || p_name == "e")
    p = pressure_function_t (pressure_function_type_t::exponential);
  else
    return error_io::incorrect_arguments;

  const std::string mode_name = argv[7];
  if (mode_name == "debug" || mode_name == "DEBUG" || mode_name == "d")
    DEBUG = true;
  else
    DEBUG = false;

  if (!is_correct ())
    return error_io::incorrect_arguments;

  functions = initial_conditions_t (programm_task);
  functions.m_mu = mu;
  functions.m_p = p;

  return error_io::success;
}
