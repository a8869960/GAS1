#include "thomas_method.h"

#define MIN_DIVISION 1e-50

thomas_method_t::thomas_method_t (const unsigned int size, const bool debug): m_size (size), DEBUG(debug)
{
  allocate_memory ();
}

thomas_method_t::thomas_method_t (const tridiagonal_matrix_t &A, const bool debug): thomas_method_t (A.get_size (), debug) {}

thomas_method_t::~thomas_method_t ()
{
  if (m_alpha) delete [] m_alpha;
  if (m_beta) delete [] m_beta;
}

error_io thomas_method_t::solve_system (const tridiagonal_matrix_t &A, double *f)
{
  const error_io status = compute_alpha_and_beta (A, f);
  if (status == error_io::division_by_zero)
    return status;
  if (DEBUG) print_error_io (status);
  assert (A.get_size () == m_size);
  double *a = A.m_a;
  double *c = A.m_c;

  const double denom = a[m_size - 1] + c[m_size - 1] * m_alpha[m_size - 1];
  if (fabs (denom) < MIN_DIVISION)
    return error_io::division_by_zero;

  A.m_x[m_size - 1] = (f[m_size - 1] - c[m_size - 1] * m_beta[m_size - 1]) / denom;
  for (unsigned int i = m_size - 1; i > 0; i--)
    {
      A.m_x[i - 1] = m_alpha[i] * A.m_x[i] + m_beta[i];
      if (fabs (A.m_x[i - 1]) < MIN_DIVISION) A.m_x[i - 1] = 0.;
    }
  return error_io::success;
}

error_io thomas_method_t::allocate_memory ()
{
  assert (m_size > 0);

  m_alpha = new double [m_size];
  m_beta = new double [m_size];

  memset (m_alpha, 0., sizeof (double) * m_size);
  memset (m_beta, 0., sizeof (double) * m_size);

  if (m_alpha == nullptr || m_beta == nullptr)
    return error_io::error_memory;

  return error_io::success;
}

error_io thomas_method_t::check_result (const tridiagonal_matrix_t &A, const bool print_error) const // compute Ax and check error
{
  const unsigned int size = A.get_size ();

  // Compute rhs0 := Ax
  double *rhs0 = new double[size];
  for (unsigned int row_num = 0; row_num < size; row_num++)
    {
      if (row_num == 0)
        rhs0[row_num] = A.m_a[row_num] * A.m_x[row_num] + A.m_b[row_num] * A.m_x[row_num + 1];
      else if (row_num == size - 1)
        rhs0[row_num] = A.m_c[row_num] * A.m_x[row_num - 1] + A.m_a[row_num] * A.m_x[row_num];
      else
        rhs0[row_num] = A.m_c[row_num] * A.m_x[row_num - 1] + A.m_a[row_num] * A.m_x[row_num] + A.m_b[row_num] * A.m_x[row_num + 1];
    }

  // Compute error
  double error = 0., rhs_norm = 0., max_diff = 0.;
  unsigned int i_max_diff = 0;
  for (unsigned int i = 0; i < size; i++)
    {
      const double diff = fabs (A.m_rhs[i] - rhs0[i]);
      error += diff;
      if (diff > max_diff)
        {
          max_diff = diff;
          i_max_diff = i;
        }
    }
  for (unsigned int i = 0; i < size; i++)
    rhs_norm += fabs (A.m_rhs[i]);
  delete [] rhs0;

  const double relative_error = error / rhs_norm;
  if (print_error)
    printf ("COMPUTATIONAL ERROR = %10.7e, MAX_DIFF[%d] = %10.7e\n", relative_error, i_max_diff, max_diff);
  if (relative_error > 1e-10)
    {
      if (print_error)
        A.print_to_konsole ();
      return error_io::big_relative_error;
    }

  return error_io::success;
}

error_io thomas_method_t::compute_alpha_and_beta (const tridiagonal_matrix_t &A, double *f)
{
  assert (A.get_size () == m_size);
  double *a = A.m_a;
  double *b = A.m_b;
  double *c = A.m_c;

  if (fabs (a[0]) < MIN_DIVISION)
    return error_io::division_by_zero;

  m_alpha[1] = - b[0] / a[0];
  m_beta[1] = f[0] / a[0];

  if (fabs (m_alpha[1]) < MIN_DIVISION) m_alpha[1] = 0.;
  if (fabs (m_beta[1]) < MIN_DIVISION) m_beta[1] = 0.;

  for (unsigned int i = 1; i < m_size - 1; i++)
    {
      double denom = a[i] + m_alpha[i] * c[i];
      if (fabs (denom) < MIN_DIVISION)
        {
          m_alpha[i + 1] = m_beta[i + 1] = 0.;
          continue;
        }
      m_alpha[i + 1] = - b[i] / denom;
      m_beta[i + 1] = (f[i] - c[i] * m_beta[i]) / denom;

      if (fabs (m_alpha[i + 1]) < MIN_DIVISION) m_alpha[i + 1] = 0.;
      if (fabs (m_beta[i + 1]) < MIN_DIVISION) m_beta[i + 1] = 0.;
    } 

  // Check stability
  for (unsigned int i = 1; i < m_size; i++)
    if (fabs (m_alpha[i]) > 1. || fabs (m_beta[i]) > 1.)
      {
        if (DEBUG)
          {
            printf ("alpha[%d] = %10.3e beta[%d] = %10.3e\n", i, m_alpha[i], i, m_beta[i]);
          }
        return error_io::not_stable;
      }

  return error_io::success;
}
