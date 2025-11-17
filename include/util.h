// =============================================================================
// emwave-c: Small utility helpers shared across modules
// =============================================================================

#ifndef EMWAVE_UTIL_H
#define EMWAVE_UTIL_H

/* Integer clamp */
static inline int util_clamp_int(int v, int lo, int hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

/* Double clamp */
static inline double util_clamp_double(double v, double lo, double hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

#endif /* EMWAVE_UTIL_H */

