#include "matrix.h"

error_io matrix_t::allocate_matrix ()
{
  assert (m_size != 0 && m_H_matrix == nullptr && m_V_matrix == nullptr);

  m_H_matrix = new double [m_size * m_size];
  m_V_matrix = new double [m_size * m_size];

  if (m_H_matrix == nullptr || m_V_matrix == nullptr)
    return error_io::error_memory;

  return error_io::success;
}

matrix_t::matrix_t(const unsigned int size): m_size (size)
{

}

matrix_t::matrix_t (const arguments_t args)
{

}

matrix_t::~matrix_t ()
{
  if (m_H_matrix) delete [] m_H_matrix;
  if (m_V_matrix) delete [] m_V_matrix;
}

void matrix_t::fill_m_V_n (const unsigned int n)
{
  for (unsigned int i = 0; i < m_size; i++)
    {
      m_V_n[i] =
    }
}
