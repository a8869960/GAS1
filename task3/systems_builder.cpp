#include "systems_builder.h"

error_io systems_builder_t::allocate_memory ()
{
  assert (m_M > 2);

  m_H_prev_layer = new double [m_M];
  m_V_prev_layer = new double [m_M];

  if (m_H_prev_layer == nullptr || m_V_prev_layer == nullptr)
    return error_io::error_memory;

  return error_io::success;
}

systems_builder_t::systems_builder_t (
    const grid_t grid,
    const initial_conditions_t functions): m_functions(functions), m_mu(functions.m_mu)
{
  m_grid = grid;

  m_M = m_grid.get_M ();
  assert (allocate_memory () == error_io::success);
}

systems_builder_t::~systems_builder_t ()
{
  if (m_H_prev_layer) delete [] m_H_prev_layer;
  if (m_V_prev_layer) delete [] m_V_prev_layer;
}

void systems_builder_t::fill_zero_layer_arrays ()
{
  for (unsigned int i = 0; i < m_M; i++)
    {
      m_H_prev_layer[i] = m_functions.rho_0 (m_grid.m (i));
      m_V_prev_layer[i] = m_functions.u_0 (m_grid.m (i));
    }
}

void systems_builder_t::fill_H_prev_layer (const tridiagonal_matrix_t &H)
{
  for (unsigned int i = 0; i < H.get_size (); i++)
    {
      m_H_prev_layer[i] = H.m_x[i];
    }
}
void systems_builder_t::fill_V_prev_layer (const tridiagonal_matrix_t &V)
{
  for (unsigned int i = 0; i < V.get_size (); i++)
    {
      m_V_prev_layer[i] = V.m_x[i];
    }
}

void systems_builder_t::fill_H_matrix (
    tridiagonal_matrix_t &H_n,
    const unsigned int n)
{
  assert (n > 0);
  auto f0 = [this] (const unsigned int n, const unsigned int m)
  {
    return m_functions.f0 (m_grid.t (n), m_grid.m (m));
  };

  if (n == 1)
    fill_zero_layer_arrays ();
  else
    fill_H_prev_layer (H_n);

  const double tau = m_grid.get_tau ();
  const double h = m_grid.get_h ();

  // Fill equation #0 (m := 0)
  H_n.m_a[0] = 1. - tau / h * m_V_prev_layer[0];
  H_n.m_b[0] = tau / h * m_V_prev_layer[1];
  H_n.m_c[0] = 0.;
  H_n.m_rhs[0] = m_H_prev_layer[0] + tau * f0 (n - 1, 0);

  // Fill equation #1, 2,..., M - 2
  for (unsigned int m = 1; m < m_M - 1; m++)
    {
      if (m_V_prev_layer[m] < 0.)
        {
          H_n.m_a[m] = 1. - tau / h * m_V_prev_layer[m];
          H_n.m_b[m] = tau / h * m_V_prev_layer[m + 1];
          H_n.m_c[m] = 0.;
        }
      else /// if m_V_prev_layer[m] >= 0.
        {
          H_n.m_a[m] = 1. + tau / h * m_V_prev_layer[m];
          H_n.m_b[m] = 0.;
          H_n.m_c[m] = - tau / h * m_V_prev_layer[m - 1];
        }
      H_n.m_rhs[m] = tau * f0 (n - 1, m) + m_H_prev_layer[m];
    }

  //Fill equation #M - 1
   H_n.m_a[m_M - 1] = 1. + tau / h * m_V_prev_layer[m_M - 1];
   H_n.m_b[m_M - 1] = 0.;
   H_n.m_c[m_M - 1] = - tau / h * m_V_prev_layer[m_M - 2];
   H_n.m_rhs[m_M - 1] = m_H_prev_layer[m_M - 1] + tau * f0 (n - 1, m_M - 1);

//   H_n.print_to_konsole ();
}

void systems_builder_t::fill_V_matrix (
    tridiagonal_matrix_t &V_n,
    const tridiagonal_matrix_t &H_n,
    const unsigned int n)
{
  assert (n > 0);
  auto f = [this] (const unsigned int n, const unsigned int m)
  {
    return m_functions.f (m_grid.t (n), m_grid.m (m));
  };

  if (n == 1)
    fill_zero_layer_arrays ();
  else
    fill_V_prev_layer (V_n);

  const double tau = m_grid.get_tau ();
  const double h = m_grid.get_h ();

  // Fill equation #0 (m := 0) (u (t, 0) = 0)
  V_n.m_a[0] = 1.;
  V_n.m_rhs[0] = 0.;

  // Fill equation #1, 2,..., M - 2
  for (unsigned int m = 1; m < m_M - 1; m++)
    {
      if (m_V_prev_layer[m] < 0.)
        {
          V_n.m_a[m] = H_n.m_x[m] * ( 1. - tau / h * m_V_prev_layer[m]) + 2 * m_mu * tau / h / h;
          V_n.m_b[m] = tau / h * H_n.m_x[m] * m_V_prev_layer[m] - m_mu * tau / h / h;
          V_n.m_c[m] = - m_mu * tau / h / h;
        }
      else /// if m_V_prev_layer[m] >= 0.
        {
          V_n.m_a[m] = H_n.m_x[m] * ( 1. + tau / h * m_V_prev_layer[m]) + 2 * m_mu * tau / h / h;
          V_n.m_b[m] = - m_mu * tau / h / h;
          V_n.m_c[m] = - tau / h * H_n.m_x[m] * m_V_prev_layer[m] - m_mu * tau / h / h;
        }
//      if (H_n.m_x[m + 1] < 0.)
//        H_n.m_x[m + 1] = 0.;
//      if (H_n.m_x[m - 1] < 0.)
//        H_n.m_x[m - 1] = 0.;
      V_n.m_rhs[m] = tau * H_n.m_x[m] * f (n, m) + H_n.m_x[m] * m_V_prev_layer[m]
          - tau / 2. / h * (m_functions.m_p.p (H_n.m_x[m + 1]) - m_functions.m_p.p (H_n.m_x[m - 1]));
    }

  //Fill equation #M - 1 (u (t, 1) = 0)
  V_n.m_a[m_M - 1] = 1.;
  V_n.m_rhs[m_M - 1] = 0.;
}
