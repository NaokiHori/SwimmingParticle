#include <math.h>
#include <stddef.h>
#include "constant.h"
#include "memory.h"
#include "./compute_stream_function.h"
#include "./exchange_halo.h"

int compute_stream_function(
    const domain_t * const domain,
    rdft_plan_t * const rdft_plan,
    const double * const surface_velocity,
    double ** const stream_function
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double * const xf = domain->xf;
#pragma omp parallel for
  for (size_t i = 1; i <= nx + 1; i++) {
    const double x = xf[i];
    double * const s_freq = &stream_function[i][1];
    s_freq[     0] = 0.;
    s_freq[ny / 2] = 0.;
    for (size_t k = 1; k < ny / 2; k++) {
      const double pref = 0.5 * (1. - pow(x, 2.)) / pow(x, k);
      s_freq[     k] = pref * surface_velocity[     k];
      s_freq[ny - k] = pref * surface_velocity[ny - k];
    }
    rdft_exec_b(rdft_plan, s_freq);
  }
  exchange_halo(domain, stream_function);
  return 0;
}

