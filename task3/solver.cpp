#include "solver.h"

#define max_iter 5000000

error_io solver_t::allocate_memory ()
{
  const unsigned int M = m_grid.get_M ();
  const unsigned int N = m_grid.get_N ();

  assert (M > 0 && N > 0);

  m_H_result = new double [M];
  m_V_result = new double [M];

  if (m_H_result == nullptr || m_V_result == nullptr)
    return error_io::error_memory;

  return error_io::success;
}

solver_t::solver_t (const grid_t grid, const arguments_t args)
{
  m_grid = grid;
  epsilon = args.eps;
  m_functions.k = args.k;
  assert (allocate_memory () == error_io::success);

  // Set initial
  m_functions = args.functions;

  DEBUG = args.DEBUG;
}

solver_t::~solver_t ()
{
  if (m_H_result) delete [] m_H_result;
  if (m_V_result) delete [] m_V_result;
}

error_io solver_t::solve_systems ()
{
  const unsigned int M = m_grid.get_M ();

  systems_builder_t systems_builder (m_grid, m_functions);
  thomas_method_t thomas_method (M, DEBUG);

  tridiagonal_matrix_t H_matrix (M);
  tridiagonal_matrix_t V_matrix (M);

  for (unsigned int i = 0; i < M; i++)
    {
      m_H_result[i] = m_functions.rho_0 (m_grid.m(i));
      m_V_result[i] = m_functions.u_0 (m_grid.m(i));
    }

  unsigned int n;
  for (n = 1; n < max_iter; n++)
    {
      // Fill and solve H matrix
      systems_builder.fill_H_matrix (H_matrix, n);
      const error_io solve_H_status = thomas_method.solve_system (H_matrix, H_matrix.m_rhs);

      if (solve_H_status != error_io::success)
        {
          iteration = n;
          return solve_H_status;
        }

      if (DEBUG)
        {
          printf ("H_matrix, iteration %d\n", n);
          H_matrix.print_to_konsole ();
        }

      // Fill and solve V matrix
      systems_builder.fill_V_matrix (V_matrix, H_matrix, n);
      const error_io solve_V_status = thomas_method.solve_system (V_matrix, V_matrix.m_rhs);

      if (solve_V_status != error_io::success)
        {
          iteration = n;
          return solve_V_status;
        }

      if (DEBUG)
        {
          printf ("V_matrix, iteration %d\n", n);
          V_matrix.print_to_konsole ();
        }

      // Write result
      memcpy (m_V_result, V_matrix.m_x, sizeof (double) * V_matrix.get_size ());
      memcpy (m_H_result, H_matrix.m_x, sizeof (double) * H_matrix.get_size ());

      if (need_to_stop (H_matrix, V_matrix))
        {
          break;
        }
    }
  printf ("eps = %.0e k = %d T_st = %lf n_st = %d  res (n_st) = %10.3e\n", epsilon, m_functions.k, m_grid.t (n), n, res);
  return error_io::success;
}

double solver_t::compute_error_u () const
{
  double result = 0.;
  for (unsigned int i = 0; i < m_grid.get_M (); i++)
    {
      const double diff = fabs (m_functions.u (m_grid.get_T (), m_grid.m (i)) - m_V_result[i]);
      if (diff > result)
        result = diff;
    }
  return result;
}

double solver_t::compute_error_rho () const
{
  double result = 0.;
  for (unsigned int i = 0; i < m_grid.get_M (); i++)
    {
      const double diff = fabs (m_functions.rho (m_grid.get_T (), m_grid.m (i)) - m_H_result[i]);
      if (diff > result)
        result = diff;
    }
  return result;
}

bool solver_t::need_to_stop (
    const tridiagonal_matrix_t &H_n,
    const tridiagonal_matrix_t &V_n)
{
  // Compute H_tilde and V_tilde
  const unsigned int M = m_grid.get_M ();
  double s = 0.;
  for (unsigned int i = 0; i < M; i++)
    s += H_n.m_x[i];
  const double H_tilde = s / M;
  const double V_tilde = 0.;

  // Find max
  double max = 0.;
  for (unsigned int i = 0; i < M; i++)
    {
      const double H_diff = fabs (H_n.m_x[i] - H_tilde);
      const double V_diff = fabs (V_n.m_x[i] - V_tilde);

      if (H_diff > max)
        max = H_diff;
      if (V_diff > max)
        max = V_diff;
    }

  res = max;
  if (max <= epsilon)
    return true;
  return false;
}

double solver_t::get_delta_massa (const tridiagonal_matrix_t &H_n) const
{
  const unsigned int M = m_grid.get_M ();
  double H_n_m_sum = 0.;
  double H_0_m_sum = 0.;

  for (unsigned int i = 0; i < M; i++)
    {
      H_n_m_sum += H_n.m_x[i];
      H_0_m_sum += m_functions.rho_0 (m_grid.m (i));
    }
  if (fabs (H_0_m_sum) < EPS)
    return -1.;

  const double delta_massa = (H_n_m_sum - H_0_m_sum) / H_0_m_sum;
  return delta_massa;
}
