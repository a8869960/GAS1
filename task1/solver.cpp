#include "solver.h"

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
  const unsigned int N = m_grid.get_N ();

  systems_builder_t systems_builder (m_grid, m_functions);
  thomas_method_t thomas_method (M, DEBUG);

  tridiagonal_matrix_t H_matrix (M);
  tridiagonal_matrix_t V_matrix (M);

  for (unsigned int n = 1; n < N + 1; n++)
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

//      print_error_io (thomas_method.check_result (H_matrix, DEBUG));
//      print_error_io (thomas_method.check_result (V_matrix, DEBUG));
    }

  // Write result
  memcpy (m_V_result, V_matrix.m_x, sizeof (double) * V_matrix.get_size ());
  memcpy (m_H_result, H_matrix.m_x, sizeof (double) * H_matrix.get_size ());
//  V_matrix.print_result_to_konsole ();
//  H_matrix.print_result_to_konsole ();

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
