/* Copyright (C) 2023 University of Texas at El Paso

   Contributed by: Christoph Lauter

                   and the 2023 class of CS4390/5390

                   Applied Numerical Computing for Multimedia
                   Applications.

   All rights reserved.

   NO LICENSE SPECIFIED.

*/

#ifndef WIDEFLOAT_OPS_H
#define WIDEFLOAT_OPS_H

#include <stdint.h>
#include "safeinteger_ops.h"

typedef enum {
  FPCLASS_NAN      = 0,
  FPCLASS_POS_INF,
  FPCLASS_NEG_INF,
  FPCLASS_NUMBER
} widefloatclass_t;

typedef struct {
  widefloatclass_t fpclass;
  unsigned int     sign:1;
  int32_t          exponent;
  size_t           mantissa_size;
  uint64_t         *mantissa;
} widefloat_t;

#define WIDEFLOAT_OVERHEAD  ((uint64_t) 11)

void widefloat_init(widefloat_t * op, size_t n);

void widefloat_clear(widefloat_t *op);

void widefloat_set_from_scaled_integer(widefloat_t *op,
                                       int s,
                                       int64_t E,
                                       const uint64_t *m,
                                       size_t n);

void widefloat_set_from_integer(widefloat_t *op,
                                const uint64_t *m,
                                size_t n);

void widefloat_set(widefloat_t *z,
                   const widefloat_t *x);

int widefloat_get_integer(uint64_t *m,
                          size_t n,
                          const widefloat_t *x);

void widefloat_set_from_safeinteger(widefloat_t *y,
                                    const safeinteger_t *x);

int widefloat_get_safeinteger(safeinteger_t *y,
                              const widefloat_t *x);

int widefloat_cmp(const widefloat_t *x,
                  const widefloat_t *y);

void widefloat_add(widefloat_t *z,
                   const widefloat_t *x,
                   const widefloat_t *y);

void widefloat_neg(widefloat_t *z,
                   const widefloat_t *x);

void widefloat_sub(widefloat_t *z,
                   const widefloat_t *x,
                   const widefloat_t *y);

void widefloat_mul(widefloat_t *z,
                   const widefloat_t *x,
                   const widefloat_t *y);

void widefloat_div(widefloat_t *z,
                   const widefloat_t *x,
                   const widefloat_t *y);

void widefloat_rsqrt(widefloat_t *z,
                     const widefloat_t *x);

void widefloat_rcpr(widefloat_t *z,
                    const widefloat_t *x);


#endif


