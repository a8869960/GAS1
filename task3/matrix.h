#ifndef MATRIX_H
#define MATRIX_H

#include "arguments.h"
#include "grid.h"

class matrix_t
{
  unsigned int m_size;
  double *m_H_matrix;     // H ~ rho
  double *m_V_matrix;     // V ~ u
  grid_t m_grid;

  double *m_V_n;              // V (t = n, m)



public:
  error_io allocate_matrix ();

  matrix_t (const unsigned int size);
  matrix_t (const arguments_t args);
  ~matrix_t ();

private:
  void fill_m_V_n (const unsigned int n);
};

#endif // MATRIX_H
