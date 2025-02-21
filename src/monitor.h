#if !defined(MONITOR_H)
#define MONITOR_H

#include "domain.h"

extern int monitor(
    const domain_t * const domain,
    const double time,
    double ** const concentration,
    double * const surface_concentration
);

#endif // MONITOR_H
