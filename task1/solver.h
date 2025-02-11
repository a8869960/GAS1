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
public:
  double *m_H_result;       // H ~ rho
  double *m_V_result;       // V ~ u

  bool DEBUG;

public:
  error_io allocate_memory ();
  solver_t (const grid_t grid, const arguments_t args);
  ~solver_t ();

  error_io solve_systems ();
  double compute_error_u () const;
  double compute_error_rho () const;

  // Debug information
  unsigned int iteration;
};

#endif // solver_t_H
