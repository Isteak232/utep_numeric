/* Copyright (C) 2023 University of Texas at El Paso

   Contributed by: Christoph Lauter

                   and the 2023 class of CS4390/5390

                   Applied Numerical Computing for Multimedia
                   Applications.

   All rights reserved.

   NO LICENSE SPECIFIED.

*/

#ifndef INTEGER_OPS_H
#define INTEGER_OPS_H

#include <stdint.h>


void addition(uint64_t *s,
              const uint64_t *a, size_t m,
              const uint64_t *b, size_t n);

void subtraction(uint64_t *s,
                 const uint64_t *a, size_t m,
                 const uint64_t *b, size_t n);

void multiplication(uint64_t *p,
                    const uint64_t *a, size_t m,
                    const uint64_t *b, size_t n);

int division(uint64_t *q, uint64_t *r,
             const uint64_t *a, size_t m,
             const uint64_t *b, size_t n);

void divide_by_ten(uint64_t *q, unsigned int *r,
                   const uint64_t *a, size_t n);

void shift_left(uint64_t *a, size_t n, size_t k);

void shift_right(uint64_t *a, size_t n, size_t k);

int comparison(const uint64_t *a,
               const uint64_t *b,
               size_t n);

int is_zero(const uint64_t *a, size_t n);

int convert_from_decimal_string(uint64_t *a, size_t n,
                                const char *str);

void convert_to_decimal_string(char *str,
                               const uint64_t *a, size_t n);

uint64_t leading_zeros(const uint64_t *a, size_t n);

#endif


