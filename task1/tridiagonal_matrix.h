#ifndef TRIDIAGONAL_MATRIX_H
#define TRIDIAGONAL_MATRIX_H

#include "enums.h"
#include <cassert>
#include <cstring>

/// a_0 b_0  0   0   0   0  ///
/// c_1 a_1 b_1  0   0   0  ///
///  0  c_2 a_2 b_2  0   0  ///
///  0   0  c_3 a_3 b_3  0  ///
///  *   *   *   *   *   *  ///
///  0   0   0   0  c_n a_n ///

class tridiagonal_matrix_t
{
  unsigned int m_size;

public:
  double *m_a;   // diagonal
  double *m_b;   // upper diagonal
  double *m_c;   // lower diagonal

  double *m_rhs;
  double *m_x;

public:
  tridiagonal_matrix_t () = default;
  tridiagonal_matrix_t (const unsigned int size);
  ~tridiagonal_matrix_t ();

  unsigned int get_size () const { return m_size; }
  void print_to_konsole () const;
  void print_result_to_konsole () const;

private:
  error_io allocate_memory ();
};

#endif // TRIDIAGONAL_MATRIX_H
