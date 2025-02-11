#ifndef Thomas_method_t_H
#define Thomas_method_t_H

#include "enums.h"
#include "tridiagonal_matrix.h"
#include "math.h"

#include <cstring>

#define EPS 1e-16

// Solve A*y = f

class thomas_method_t
{
  unsigned int m_size;

  double *m_alpha;
  double *m_beta;

  bool DEBUG;

public:
  thomas_method_t (const unsigned int size, const bool debug = false);
  thomas_method_t (const tridiagonal_matrix_t &A, const bool debug);
  ~thomas_method_t ();

  error_io solve_system (const tridiagonal_matrix_t &A, double *f);
  error_io solve_system (const tridiagonal_matrix_t &A) { return solve_system (A, A.m_rhs); }

  error_io check_result (const tridiagonal_matrix_t &A, const bool print_error = false) const;

private:
  error_io allocate_memory ();
  error_io compute_alpha_and_beta (const tridiagonal_matrix_t &A, double *f);
};

#endif // Thomas_method_t_H
