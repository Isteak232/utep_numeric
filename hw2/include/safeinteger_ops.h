/* Copyright (C) 2023 University of Texas at El Paso

   Contributed by: Christoph Lauter

                   and the 2023 class of CS4390/5390

                   Applied Numerical Computing for Multimedia
                   Applications.

   All rights reserved.

   NO LICENSE SPECIFIED.

*/

#ifndef SAFEINTEGER_OPS_H
#define SAFEINTEGER_OPS_H

#include <stdint.h>

typedef struct {
  int              sign:1;
  size_t           size;
  uint64_t         *data;
} safeinteger_t;

void safeinteger_init(safeinteger_t *x);

void safeinteger_clear(safeinteger_t *x);

void safeinteger_set(safeinteger_t *y,
                     const safeinteger_t *x);

void safeinteger_set_from_ui(safeinteger_t *y,
                             unsigned long long int x);

void safeinteger_set_from_si(safeinteger_t *y,
                             long long int x);

void safeinteger_set_from_integer(safeinteger_t *y,
                                  uint64_t *x,
                                  size_t size);

int safeinteger_get_ui(unsigned long long int *y,
                       const safeinteger_t *x);

int safeinteger_get_si(long long int *y,
                       const safeinteger_t *x);

void safeinteger_neg(safeinteger_t *y,
                     const safeinteger_t *x);

void safeinteger_add(safeinteger_t *z,
                     const safeinteger_t *x,
                     const safeinteger_t *y);

void safeinteger_sub(safeinteger_t *z,
                     const safeinteger_t *x,
                     const safeinteger_t *y);

void safeinteger_mul(safeinteger_t *z,
                     const safeinteger_t *x,
                     const safeinteger_t *y);

int safeinteger_div(safeinteger_t *quot,
                    safeinteger_t *rema,
                    const safeinteger_t *x,
                    const safeinteger_t *y);

void safeinteger_shift_left(safeinteger_t *z,
                            const safeinteger_t *x,
                            size_t y);

void safeinteger_shift_right(safeinteger_t *z,
                             const safeinteger_t *x,
                             size_t y);

int safeinteger_cmp(const safeinteger_t *x,
                    const safeinteger_t *y);

int safeinteger_from_string(safeinteger_t *x,
                            const char *str);

char *safeinteger_to_string(const safeinteger_t *x);

void safeinteger_to_array(int *sign,
                          uint32_t **abs,
                          size_t *size,
                          const safeinteger_t *x);

void safeinteger_from_array(safeinteger_t *x,
                            int sign,
                            uint32_t *abs,
                            size_t size);

#endif
