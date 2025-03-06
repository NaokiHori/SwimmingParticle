#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include "./compute_stream_function.h"
#include "./exchange_halo.h"

static int get_f(
    const size_t k,
    const double x,
    double * const f
) {
  // general solution
  if (1 == k) {
    f[0] = pow(x, -1.);
    f[1] = x;
    f[2] = x * log(x);
    f[3] = pow(x, 3.);
  } else {
    const double dk = (double)k;
    f[0] = pow(x, - dk);
    f[1] = pow(x, + dk);
    f[2] = pow(x, - dk + 2.);
    f[3] = pow(x, + dk + 2.);
  }
  return 0;
}

static int get_g(
    const size_t k,
    const double x,
    double * const g
) {
  // derivative of general solution
  if (1 == k) {
    g[0] = - pow(x, -2.);
    g[1] = 1.;
    g[2] = log(x) + 1.;
    g[3] = 3. * pow(x, 2.);
  } else {
    const double dk = (double)k;
    g[0] = - dk * pow(x, - dk - 1.);
    g[1] = + dk * pow(x, + dk - 1.);
    g[2] = (- dk + 2.) * pow(x, - dk + 1.);
    g[3] = (+ dk + 2.) * pow(x, + dk + 1.);
  }
  return 0;
}

static int solve_linear_systems(
    double * const a,
    double * const real,
    double * const imag
) {
  // solve 4-by-4 linear systems
  const size_t nitems = 4;
  // forward sweep
  for (size_t i = 0; i < nitems; i++) {
    // normalize i-th row
    const double f = 1. / a[i * nitems + i];
    for (size_t j = i; j < nitems; j++) {
      a[i * nitems + j] *= f;
    }
    real[i] *= f;
    imag[i] *= f;
    // eliminate lower-triangular part
    for (size_t ii = i + 1; ii < nitems; ii++) {
      const double f = a[ii * nitems + i];
      for (size_t j = i + 1; j < nitems; j++) {
        a[ii * nitems + j] -= f * a[i * nitems + j];
      }
      real[ii] -= f * real[i];
      imag[ii] -= f * imag[i];
    }
  }
  // backward substitution
  for (size_t i = nitems - 1; ; i--) {
    for (size_t j = i + 1; j < nitems; j++) {
      real[i] -= a[i * nitems + j] * real[j];
      imag[i] -= a[i * nitems + j] * imag[j];
    }
    if (0 == i) {
      break;
    }
  }
  return 0;
}

static double dot(
    const double a[4],
    const double b[4]
) {
  // compute inner product
  double result = 0.;
  for (size_t n = 0; n < 4; n++) {
    result += a[n] * b[n];
  }
  return result;
}

int compute_stream_function(
    const domain_t * const domain,
    rdft_plan_t * const rdft_plan,
    const double * const surface_velocity,
    double ** const stream_function
) {
  // impose no-slip and impermable condition rigorously
  //   on the outer wall or not
  const bool use_rigorous = false;
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double * const xf = domain->xf;
  // from surface velocity (and other boundary conditions),
  //   compute stream function in the azimuthal frequency space
  //   for each radial position
  // k = 0 and ny / 2: the surface velocity is zero and thus zero
  for (size_t i = 1; i <= nx + 1; i++) {
    stream_function[i][         1] = 0.;
    stream_function[i][ny / 2 + 1] = 0.;
  }
  // other azimuthal wavenumbers
  if (use_rigorous) {
#pragma omp parallel for
    for (size_t k = 1; k < ny / 2; k++) {
      // solve a 4-by-4 linear system to find coefficients
      //   at this wavenumber
      double matrix[16] = {0.};
      get_f(k, xf[     1], matrix +  0);
      get_f(k, xf[nx + 1], matrix +  4);
      get_g(k, xf[     1], matrix +  8);
      get_g(k, xf[nx + 1], matrix + 12);
      // real / imaginary parts are handled separately
      double coefs_real[4] = {
        0.,
        0.,
        - surface_velocity[k],
        0.,
      };
      double coefs_imag[4] = {
        0.,
        0.,
        + surface_velocity[ny - k],
        0.,
      };
      solve_linear_systems(matrix, coefs_real, coefs_imag);
      for (size_t i = 1; i <= nx + 1; i++) {
        const double x = xf[i];
        double buf[4] = {0.};
        get_f(k, x, buf);
        double * const real = &stream_function[i][     k + 1];
        double * const imag = &stream_function[i][ny - k + 1];
        *real = + dot(coefs_real, buf);
        *imag = - dot(coefs_imag, buf);
      }
    }
  } else {
    // Hu et al., JCP, 2019
#pragma omp parallel for
    for (size_t k = 1; k < ny / 2; k++) {
      for (size_t i = 1; i <= nx + 1; i++) {
        const double x = xf[i];
        const double pref = 0.5 * (1. - pow(x, 2.)) / pow(x, k);
        double * const real = &stream_function[i][     k + 1];
        double * const imag = &stream_function[i][ny - k + 1];
        *real = pref * surface_velocity[     k];
        *imag = pref * surface_velocity[ny - k];
      }
    }
  }
  // transform to physical space
#pragma omp parallel for
  for (size_t i = 1; i <= nx + 1; i++) {
    rdft_exec_b(rdft_plan, &stream_function[i][1]);
  }
  exchange_halo(domain, stream_function);
  return 0;
}

