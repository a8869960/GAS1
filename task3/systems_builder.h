#ifndef SYSTEMS_BUILDER_H
#define SYSTEMS_BUILDER_H

#include "grid.h"
#include "tridiagonal_matrix.h"
#include "initial_conditions.h"

class systems_builder_t
{
  // Initial conditions
  initial_conditions_t m_functions;

  double *m_H_prev_layer;           // m_H_prev_layer[m] := H^{n - 1}_m
  double *m_V_prev_layer;           // m_V_prev_layer[m] := V^{n - 1}_m

  unsigned int m_M;                 // size of m_H_0_m

  double m_mu;

  grid_t m_grid;
public:
  systems_builder_t () = default;
  systems_builder_t (const grid_t grid, const initial_conditions_t functions);
  ~systems_builder_t ();

  void fill_H_matrix (
      tridiagonal_matrix_t &H_n,
      const unsigned int n);  // fill H matrix on n step
  void fill_V_matrix (
      tridiagonal_matrix_t &V_n,
      const tridiagonal_matrix_t &H_n,
      const unsigned int n);  // fill /// a_0 b_0 matrix on n step


private:
  error_io allocate_memory ();
  void fill_zero_layer_arrays ();    /// fill m_H_prev_layer and m_V_prev_layer
  void fill_H_prev_layer (const tridiagonal_matrix_t &H);
  void fill_V_prev_layer (const tridiagonal_matrix_t &V);
};

#endif // SYSTEMS_BUILDER_H
