#if !defined(COMPUTE_VELOCITY_H)
#define COMPUTE_VELOCITY_H

#include "domain.h"
#include "rdft.h"

extern int compute_stream_function(
    const domain_t * const domain,
    rdft_plan_t * const rdft_plan,
    const double * const surface_concentration,
    double ** const stream_function
);

#endif // COMPUTE_VELOCITY_H
