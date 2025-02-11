#include "tridiagonal_matrix.h"

tridiagonal_matrix_t::tridiagonal_matrix_t (const unsigned int size): m_size(size)
{
  assert (allocate_memory () == error_io::success);
}

tridiagonal_matrix_t::~tridiagonal_matrix_t ()
{
  if (m_a) delete [] m_a;
  if (m_b) delete [] m_b;
  if (m_c) delete [] m_c;
  if (m_rhs) delete [] m_rhs;
  if (m_x) delete [] m_x;
}

void tridiagonal_matrix_t::print_to_konsole () const
{
  if (m_size > 15)
    return;

  printf ("\nSYSTEM:\n");
  // Print #0 row
  for (unsigned int i = 0; i < m_size; i++)
    {
      if (i == 0)
        printf (" %10.3e", m_a[0]);
      else if (i == 1)
        printf (" %10.3e", m_b[0]);
      else
        printf (" %10.3e", 0.);
    }
  printf (" | %10.3e\n", m_rhs[0]);

  // Print #1 ... m_size - 2 row
  for (unsigned int row_num = 1; row_num < m_size - 1; row_num++)
    {
      const unsigned int shift = row_num - 1;
      for (unsigned int i = 0; i < shift; i++)
        printf (" %10.3e", 0.);
      printf (" %10.3e %10.3e %10.3e", m_c[row_num], m_a[row_num], m_b[row_num]);
      for (unsigned int i = shift + 3; i < m_size; i++)
        printf (" %10.3e", 0.);
      printf (" | %10.3e\n", m_rhs[row_num]);
    }

  // Print #m_size - 1 row
  for (unsigned int i = 0; i < m_size; i++)
    {
      if (i == m_size - 1)
        printf (" %10.3e", m_a[m_size - 1]);
      else if (i == m_size - 2)
        printf (" %10.3e", m_c[m_size - 1]);
      else
        printf (" %10.3e", 0.);
    }
  printf (" | %10.3e\n", m_rhs[m_size - 1]);
  print_result_to_konsole ();
}

void tridiagonal_matrix_t::print_result_to_konsole () const
{
  printf ("RESULT:\n");
  for (unsigned int i = 0; i < m_size; i++)
    printf (" %10.3e", m_x[i]);
  printf ("\n");
}

error_io tridiagonal_matrix_t::allocate_memory ()
{
  assert (m_size > 0);
//  assert (m_a == nullptr && m_b == nullptr && m_c == nullptr && m_rhs == nullptr && m_x == nullptr);

  m_a = new double [m_size];
  m_b = new double [m_size];
  m_c = new double [m_size];

  m_rhs = new double [m_size];
  m_x = new double [m_size];

  if (m_a == nullptr || m_b == nullptr || m_c == nullptr || m_rhs == nullptr || m_x == nullptr)
    return error_io::error_memory;

  memset (m_a, 0., sizeof (double) * m_size);
  memset (m_b, 0., sizeof (double) * m_size);
  memset (m_c, 0., sizeof (double) * m_size);
  memset (m_rhs, 0., sizeof (double) * m_size);
  memset (m_x, 0., sizeof (double) * m_size);

  return error_io::success;
}
