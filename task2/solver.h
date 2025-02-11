#ifndef solver_t_H
#define solver_t_H

#include "tridiagonal_matrix.h"
#include "systems_builder.h"
#include "thomas_method.h"
#include "arguments.h"
#include "initial_conditions.h"

#include <cstring>

class solver_t
{
  // Initial conditions
  initial_conditions_t m_functions;

  grid_t m_grid;
  double *m_H_result;       // H ~ rho
  double *m_V_result;       // V ~ u

  bool DEBUG;

  unsigned int n_st;
  double epsilon;

public:
  error_io allocate_memory ();
  solver_t (const grid_t grid, const arguments_t args);
  ~solver_t ();

  error_io solve_systems ();
  double compute_error_u () const;
  double compute_error_rho () const;

  // Debug information
  unsigned int iteration;
private:
  bool need_to_stop (
      const tridiagonal_matrix_t &H_n,
      const tridiagonal_matrix_t &V_n);
  double get_delta_massa (const tridiagonal_matrix_t &H_n) const;

public:
  double T;
  double res;
};

#endif // solver_t_H
