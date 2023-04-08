
/* Copyright (C) 2023 University of Texas at El Paso

   Contributed by: Christoph Lauter

                   and the 2023 class of CS4390/5390

                   Applied Numerical Computing for Multimedia
                   Applications.

   All rights reserved.

   NO LICENSE SPECIFIED.

*/

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include "integer_ops.h"
#include "safeinteger_ops.h"
#include "widefloat_ops.h"

/* Helper functions */

/* Tries to multiply the two size_t arguments a and b.

   If the product holds on a size_t variable, sets the
   variable pointed to by c to that product and returns a
   non-zero value.

   Otherwise, does not touch the variable pointed to by c and
   returns zero.

   This implementation is kind of naive as it uses a division.
   If performance is an issue, try to speed it up by avoiding
   the division while making sure that it still does the right
   thing (which is hard to prove).

*/
static inline int __try_size_t_multiply(size_t *c, size_t a, size_t b) {
  size_t t, r, q, M;

  /* If any of the arguments a and b is zero, everthing works just fine. */
  if ((a == ((size_t) 0)) ||
      (b == ((size_t) 0))) {
    *c = a * b;
    return 1;
  }

  /* If both a and b are less than 2^(k/2), where k is the bitwith of
     a size_t, a regular multiplication is enough.
  */
  M = ((size_t) 1) << (((size_t) 4) * sizeof(size_t));
  if ((a < M) && (b < M)) {
    *c = a * b;
    return 1;
  }

  /* Here, neither a nor b is zero.

     We perform the multiplication, which may overflow, i.e. present
     some modulo-behavior.

  */
  t = a * b;

  /* Perform Euclidian division on t by a:

     t = a * q + r

     As we are sure that a is non-zero, we are sure
     that we will not divide by zero.

  */
  q = t / a;
  r = t % a;

  /* If the rest r is non-zero, the multiplication overflowed. */
  if (r != ((size_t) 0)) return 0;

  /* Here the rest r is zero, so we are sure that t = a * q.

     If q is different from b, the multiplication overflowed.
     Otherwise we are sure that t = a * b.

  */
  if (q != b) return 0;
  *c = t;
  return 1;
}

static inline void __m_memset(void *s, int c, size_t m, size_t n) {
  size_t p, i, min_m_n, max_m_n;
  void *curr;

  /* Easy case for 99.9999% of all cases */
  if (__try_size_t_multiply(&p, m, n)) {
    memset(s, c, p);
    return;
  }

  /* Overflow case */
  if (m < n) {
    min_m_n = m;
    max_m_n = n;
  } else {
    min_m_n = n;
    max_m_n = m;
  }
  for (i=0,curr=s;
       i<min_m_n;
       i++,curr=(void *) (((char *) curr) + max_m_n)) {
    memset(curr, c, max_m_n);
  }
}

static inline void __m_memcpy(void *dst, const void *src, size_t m, size_t n) {
  size_t p, i, min_m_n, max_m_n;
  void *curr_dst;
  const void *curr_src;

  /* Easy case for 99.9999% of all cases */
  if (__try_size_t_multiply(&p, m, n)) {
    memcpy(dst, src, p);
    return;
  }

  /* Overflow case */
  if (m < n) {
    min_m_n = m;
    max_m_n = n;
  } else {
    min_m_n = n;
    max_m_n = m;
  }
  for (i=0,curr_dst=dst,curr_src=src;
       i<min_m_n;
       i++,
         curr_dst=(void *) (((char *) curr_dst) + max_m_n),
         curr_src=(const void *) (((const char *) curr_src) + max_m_n)) {
    memcpy(curr_dst, curr_src, max_m_n);
  }
}

static inline void *__alloc_mem(size_t nmemb, size_t size) {
  void *ptr;

  ptr = calloc(nmemb, size);
  if (ptr == NULL) {
    fprintf(stderr, "Cannot allocate memory: %s\n", strerror(errno));
    exit(1);
  }

  return ptr;
}

static inline void __free_mem(void *ptr) {
  free(ptr);
}


/* Functions to (de-)allocate memory for wide floating-point
   numbers
*/

/* Initializes the widefloat_t op to NaN, allocating memory for a
   mantissa of 64 * n - WIDEFLOAT_OVERHEAD bits of precision.

   Does nothing if n is zero.

*/
void widefloat_init(widefloat_t * op, size_t n) {
  if (n == ((size_t) 0)) return;
  op->fpclass = FPCLASS_NAN;
  op->sign = 0;
  op->exponent = (int32_t) 0;
  op->mantissa_size = n;
  op->mantissa = __alloc_mem(n, sizeof(*(op->mantissa)));
}


/* Deallocates the memory in the mantissa of widefloat_t op */
void widefloat_clear(widefloat_t *op) {
  op->fpclass = FPCLASS_NAN;
  op->sign = 0;
  op->exponent = (int32_t) 0;
  op->mantissa_size = (size_t) 0;
  __free_mem(op->mantissa);
  op->mantissa = NULL;
}



/* General rounding function

   Sets the floating-point number op to the value equal to or in
   magnitude just below

   (-1)^s * 2^E * m

   where

   m is an integer with n digits.

   The floating-point number op must be initialized.

   Does nothing if op is clearly not initialized.

   If n is zero, sets op to zero.

   If E is too great, sets op to infinity.
   If E is too small, sets op to zero.

*/
void widefloat_set_from_scaled_integer(widefloat_t *op,
                                       int s,
                                       int64_t E,
                                       const uint64_t *m,
                                       size_t n) {
  uint64_t *t;
  int64_t EE;
  size_t q;
  uint64_t lzc;
  size_t sigma;
  int32_t expo;

  /* Check special cases */
  if (op->mantissa_size == ((size_t) 0)) return;
  if (op->mantissa == NULL) return;
  if (n == ((size_t) 0)) {
    op->fpclass = FPCLASS_NUMBER;
    op->sign = !!s;
    op->exponent = (int32_t) 0;
    __m_memset(op->mantissa, 0, op->mantissa_size, sizeof(*(op->mantissa)));
    return;
  }
  if (is_zero(m, n)) {
    op->fpclass = FPCLASS_NUMBER;
    op->sign = !!s;
    op->exponent = (int32_t) 0;
    __m_memset(op->mantissa, 0, op->mantissa_size, sizeof(*(op->mantissa)));
    return;
  }

  /* Here m has at least one digit and is not zero. The floating-point
     number is initialized.

     We copy m into a temporary t with size q = n + 1, so that we can
     normalize it.

  */
  q = n + ((size_t) 1);
  t = __alloc_mem(q, sizeof(*t));
  t[q - ((size_t) 1)] = (uint64_t) 0;
  __m_memcpy(t, m, n, sizeof(*t));

  /* We get a leading-zero count on t */
  lzc = leading_zeros(t, q);

  /* We shift t to the left by lzc - WIDEFLOAT_OVERHEAD */
  sigma = ((size_t) lzc) - ((size_t) WIDEFLOAT_OVERHEAD);
  shift_left(t, q, sigma);

  /* We adapt EE so that 2^EE * t = 2^E * m */
  EE = E - ((int64_t) sigma);

  /* We adapt EE to reflect a mantissa between 1 and 2 */
  EE += ((int64_t) (q << 6)) -
    ((int64_t) WIDEFLOAT_OVERHEAD) -
    ((int64_t) 1);

  /* Now we check if EE holds on a 32 bit signed integer.

     If EE is greater than the greatest 32 bit signed integer,
     we produce infinity.

     If EE is less than the smallest 32 bit signed integer,
     we produce zero.

  */
  if (EE > ((int64_t) ((((uint64_t) 1) << 31) - ((uint64_t) 1)))) {
    /* Produce infinity */
    if (s) {
      op->fpclass = FPCLASS_NEG_INF;
    } else {
      op->fpclass = FPCLASS_POS_INF;
    }
    op->sign = !!s;
    op->exponent = (int64_t) 0;
    __m_memset(op->mantissa, 0, op->mantissa_size, sizeof(*(op->mantissa)));
    __free_mem(t);
    return;
  }
  if (EE < ((-((int64_t) ((((uint64_t) 1) << 31) - ((uint64_t) 1)))) - ((int64_t) 1))) {
    /* Produce zero */
    op->fpclass = FPCLASS_NUMBER;
    op->sign = !!s;
    op->exponent = (int32_t) 0;
    __m_memset(op->mantissa, 0, op->mantissa_size, sizeof(*(op->mantissa)));
    __free_mem(t);
    return;
  }

  /* Now we know that the exponent holds on a 32bit signed integer */
  expo = (int32_t) EE;

  /* Now we can store the sign, the exponent and the mantissa in
     the floating-point number.

     As we do round-to-zero, rounding is just a truncation.

  */
  op->fpclass = FPCLASS_NUMBER;
  op->sign = !!s;
  op->exponent = expo;
  if (op->mantissa_size >= q) {
    /* The mantissa is longer than the temporary t */
    __m_memset(op->mantissa, 0, op->mantissa_size, sizeof(*(op->mantissa)));
    __m_memcpy(&(op->mantissa[op->mantissa_size - q]), t,
               q, sizeof(*(op->mantissa)));
  } else {
    /* The mantissa is shorter than the temporary t */
    __m_memcpy(op->mantissa, &t[q - op->mantissa_size],
               op->mantissa_size, sizeof(*(op->mantissa)));
  }
  __free_mem(t);
}


/* Set a floating-point number from an integer

   Sets the floating-point number op to the value equal to or in
   magnitude just below

   m

   where

   m is an integer with n digits.

   The floating-point number op must be initialized.

   Does nothing if op is clearly not initialized.

   If n is zero, sets op to zero.

*/
void widefloat_set_from_integer(widefloat_t *op,
                                const uint64_t *m,
                                size_t n) {
  widefloat_set_from_scaled_integer(op, 0, (int64_t) 0, m, n);
}


/* Floating-point multiplication

   z ~= x * y

   Special cases:

   x = NaN => return NaN

   y = NaN => return NaN

   x = +/- Inf
       if y = 0 => return NaN
       otherwise => return +/- Inf resp. -/+ Inf

   y = +/- Inf
       if x = 0 => return NaN
       otherwise => return +/- Inf

   otherwise we compute

*/
void widefloat_mul(widefloat_t *z,
                   const widefloat_t *x,
                   const widefloat_t *y) {
  int signx, signy;
  uint64_t *R;

  /* Check for NaN in input */
  if ((x->fpclass == FPCLASS_NAN) ||
      (y->fpclass == FPCLASS_NAN)) {
    z->fpclass = FPCLASS_NAN;
    z->sign = 0;
    z->exponent = 0;
    __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
    return;
  }

  /* Check if x = Inf */
  if ((x->fpclass == FPCLASS_POS_INF) ||
      (x->fpclass == FPCLASS_NEG_INF)) {
    /* Check if y is zero */
    if ((y->fpclass == FPCLASS_NUMBER) &&
        is_zero(y->mantissa, y->mantissa_size)) {
      /* x = +/- Inf, y = +/- 0

         Return NaN.

      */
      z->fpclass = FPCLASS_NAN;
      z->sign = 0;
      z->exponent = 0;
      __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
      return;
    }
    signy = y->sign;
    if (y->fpclass == FPCLASS_POS_INF) signy = 0;
    if (y->fpclass == FPCLASS_NEG_INF) signy = 1;
    if (((x->fpclass == FPCLASS_POS_INF) && (!!signy == 0)) ||
        ((x->fpclass == FPCLASS_NEG_INF) && (!!signy == 1))) {
      z->fpclass = FPCLASS_POS_INF;
      z->sign = 0;
    } else {
      z->fpclass = FPCLASS_NEG_INF;
      z->sign = 1;
    }
    z->exponent = 0;
    __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
    return;
  }

  /* Check if y = Inf */
  if ((y->fpclass == FPCLASS_POS_INF) ||
      (y->fpclass == FPCLASS_NEG_INF)) {
    /* Check if x is zero */
    if ((x->fpclass == FPCLASS_NUMBER) &&
        is_zero(x->mantissa, x->mantissa_size)) {
      /* y = +/- Inf, x = +/- 0

         Return NaN.

      */
      z->fpclass = FPCLASS_NAN;
      z->sign = 0;
      z->exponent = 0;
      __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
      return;
    }
    signx = x->sign;
    if (x->fpclass == FPCLASS_POS_INF) signx = 0;
    if (x->fpclass == FPCLASS_NEG_INF) signx = 1;
    if (((y->fpclass == FPCLASS_POS_INF) && (!!signx == 0)) ||
        ((y->fpclass == FPCLASS_NEG_INF) && (!!signx == 1))) {
      z->fpclass = FPCLASS_POS_INF;
      z->sign = 0;
    } else {
      z->fpclass = FPCLASS_NEG_INF;
      z->sign = 1;
    }
    z->exponent = 0;
    __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
    return;
  }

  /* Here x and y are numbers.

     We multiply their mantissas.

     We need to allocate enough space first.

  */
  R = __alloc_mem(x->mantissa_size + y->mantissa_size, sizeof(*R));
  multiplication(R,
                 x->mantissa, x->mantissa_size,
                 y->mantissa, y->mantissa_size);
  widefloat_set_from_scaled_integer(z,
                                    !!(x->sign ^ y->sign),
                                    (((int64_t) x->exponent) +
                                     ((int64_t) y->exponent) +
                                     ((int64_t) 2) +
                                     (((int64_t) 2) * ((int64_t) WIDEFLOAT_OVERHEAD)) -
                                     (((int64_t) 64) * (((int64_t) x->mantissa_size) +
                                                        ((int64_t) y->mantissa_size)))),
                                    R,
                                    x->mantissa_size + y->mantissa_size);
  __free_mem(R);
}


/* Floating-point addition

   z ~= x + y

   Special cases:

   x = NaN => return NaN

   y = NaN => return NaN

   x = 0 => return y

   y = 0 => return x

   x = +/- Inf
       if y = -/+ Inf => NaN
       otherwise => +/- Inf

   y = +/- Inf
       if x = -/+ Inf => NaN
       otherwise => +/- Inf

   otherwise we compute

*/
void widefloat_add(widefloat_t *z,
                   const widefloat_t *x,
                   const widefloat_t *y) {
  uint64_t *R;
  uint64_t *M;
  uint64_t *N;
  int64_t E;
  int64_t F;
  int s;
  int t;
  uint64_t p;
  uint64_t q;
  uint64_t r;
  size_t pr;
  size_t qr;
  uint64_t *SM;
  uint64_t *SN;
  size_t SMs;
  size_t SNs;
  size_t Ks;
  uint64_t *K;
  int64_t sigma;
  int64_t G;
  uint64_t *SMM;
  uint64_t *SNN;

  /* Check for NaN in input */
  if ((x->fpclass == FPCLASS_NAN) ||
      (y->fpclass == FPCLASS_NAN)) {
    z->fpclass = FPCLASS_NAN;
    z->sign = 0;
    z->exponent = 0;
    __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
    return;
  }

  /* Check if x is 0 */
  if ((x->fpclass == FPCLASS_NUMBER) &&
      is_zero(x->mantissa, x->mantissa_size)) {
    /* Return y */
    if (y->fpclass == FPCLASS_NUMBER) {
      /* y is a number. The precision of z can be
         different. So we can have a rounding or appending
         of zeros.
      */
      R = __alloc_mem(y->mantissa_size, sizeof(*R));
      __m_memcpy(R, y->mantissa, y->mantissa_size, sizeof(*R));
      widefloat_set_from_scaled_integer(z,
                                        y->sign,
                                        (((int64_t) y->exponent) +
                                         ((int64_t) 1) -
                                         (((int64_t) y->mantissa_size) *
                                          ((int64_t) 64)) +
                                         ((int64_t) WIDEFLOAT_OVERHEAD)),
                                        R,
                                        y->mantissa_size);
      __free_mem(R);
    } else {
      /* y is Inf (or NaN) */
      z->fpclass = y->fpclass;
      z->sign = 0;
      z->exponent = 0;
      __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
    }
    return;
  }

  /* Check if y is 0 */
  if ((y->fpclass == FPCLASS_NUMBER) &&
      is_zero(y->mantissa, y->mantissa_size)) {
    /* Return x */
    if (x->fpclass == FPCLASS_NUMBER) {
      /* x is a number. The precision of z can be
         different. So we can have a rounding or appending
         of zeros.
      */
      R = __alloc_mem(x->mantissa_size, sizeof(*R));
      __m_memcpy(R, x->mantissa, x->mantissa_size, sizeof(*R));
      widefloat_set_from_scaled_integer(z,
                                        x->sign,
                                        (((int64_t) x->exponent) +
                                         ((int64_t) 1) -
                                         (((int64_t) x->mantissa_size) *
                                          ((int64_t) 64)) +
                                         ((int64_t) WIDEFLOAT_OVERHEAD)),
                                        R,
                                        x->mantissa_size);
      __free_mem(R);
    } else {
      /* x is Inf (or NaN) */
      z->fpclass = x->fpclass;
      z->sign = 0;
      z->exponent = 0;
      __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
    }
    return;
  }

  /* Check if x = +/- Inf */
  if ((x->fpclass == FPCLASS_POS_INF) ||
      (x->fpclass == FPCLASS_NEG_INF)) {
    /* Check if y is the other infinity */
    if (((x->fpclass == FPCLASS_POS_INF) &&
         (y->fpclass == FPCLASS_NEG_INF)) ||
        ((x->fpclass == FPCLASS_NEG_INF) &&
         (y->fpclass == FPCLASS_POS_INF))) {
      /* Return NaN */
      z->fpclass = FPCLASS_NAN;
      z->sign = 0;
      z->exponent = 0;
      __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
      return;
    }
    /* Return x = +/- Inf */
    z->fpclass = x->fpclass;
    z->sign = 0;
    z->exponent = 0;
    __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
    return;
  }

  /* Check if y = +/- Inf */
  if ((y->fpclass == FPCLASS_POS_INF) ||
      (y->fpclass == FPCLASS_NEG_INF)) {
    /* Check if x is the other infinity */
    if (((y->fpclass == FPCLASS_POS_INF) &&
         (x->fpclass == FPCLASS_NEG_INF)) ||
        ((y->fpclass == FPCLASS_NEG_INF) &&
         (x->fpclass == FPCLASS_POS_INF))) {
      /* Return NaN */
      z->fpclass = FPCLASS_NAN;
      z->sign = 0;
      z->exponent = 0;
      __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
      return;
    }
    /* Return y = +/- Inf */
    z->fpclass = y->fpclass;
    z->sign = 0;
    z->exponent = 0;
    __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
    return;
  }

  /* Here x and y are numbers and none of both is zero

     Order x and y so that x has the greater exponent.

     If x and y have the same exponent, order them
     in such a way that the mantissa of x is greater
     than that of y.

  */
  r = ((uint64_t) z->mantissa_size) * ((uint64_t) 64)
    - ((uint64_t) WIDEFLOAT_OVERHEAD);
  if (x->exponent >= y->exponent) {
    s = x->sign;
    t = y->sign;
    E = (int64_t) x->exponent;
    F = (int64_t) y->exponent;
    p = ((uint64_t) x->mantissa_size) * ((uint64_t) 64)
      - ((uint64_t) WIDEFLOAT_OVERHEAD);
    q = ((uint64_t) y->mantissa_size) * ((uint64_t) 64)
      - ((uint64_t) WIDEFLOAT_OVERHEAD);
    M = x->mantissa;
    N = y->mantissa;
    pr = x->mantissa_size;
    qr = y->mantissa_size;
  } else {
    s = y->sign;
    t = x->sign;
    E = (int64_t) y->exponent;
    F = (int64_t) x->exponent;
    p = ((uint64_t) y->mantissa_size) * ((uint64_t) 64)
      - ((uint64_t) WIDEFLOAT_OVERHEAD);
    q = ((uint64_t) x->mantissa_size) * ((uint64_t) 64)
      - ((uint64_t) WIDEFLOAT_OVERHEAD);
    M = y->mantissa;
    N = x->mantissa;
    pr = y->mantissa_size;
    qr = x->mantissa_size;
  }

  /* We need to compute

     (-1)^s * 2^E * 2^(-p+1) * M + (-1)^t * 2^F * 2^(-q+1) * N

     where E >= F.

  */
  if (F < (E - ((int64_t) 5) - ((int64_t) r))) {
    /* (-1)^t * 2^F * 2^(-q+1) * N

       is so small with respect to

       (-1)^s * 2^E * 2^(-p+1) * M

       so that we can just return

       (-1)^s * 2^E * 2^(-p+1) * M

    */
    R = __alloc_mem(pr, sizeof(*R));
    __m_memcpy(R, M, pr, sizeof(*R));
    widefloat_set_from_scaled_integer(z,
                                      s,
                                      (E
                                       + ((int64_t) 1)
                                       - ((int64_t) p)),
                                      R,
                                      pr);
    __free_mem(R);
    return;
  }

  /* E - F <= r + 5

     Let be

     sigma = E - F - p + q


  */
  sigma = E - F - ((int64_t) p) + ((int64_t) q);
  if (sigma >= ((int64_t) 0)) {
    G = F - ((int64_t) q) + ((int64_t) 1);
    /* Shift M sigma to the left,
       copy N

       Allocate SM as size of M + sigma / 64 + 1
       Allocate SN as size of N

    */
    SMs = pr
      + ((size_t) (((uint64_t) sigma) / ((uint64_t) 64)))
      + ((size_t) 1);
    SNs = qr;
    SM = __alloc_mem(SMs, sizeof(*SM));
    SN = __alloc_mem(SNs, sizeof(*SN));
    __m_memset(SM, 0,
               SMs,
               sizeof(*SM));
    __m_memcpy(SM, M, pr, sizeof(*SM));
    __m_memcpy(SN, N, qr, sizeof(*SN));
    shift_left(SM, SMs, (size_t) sigma);
  } else {
    G = F - ((int64_t) q) + ((int64_t) 1) - ((int64_t) sigma);
    /* Shift N -sigma to the left and copy M

       Allocate SM as size of M
       Allocate SN as size of N + (-sigma) / 64 + 1

    */
    SMs = pr;
    SNs = qr
      + ((size_t) (((uint64_t) (-sigma)) / ((uint64_t) 64)))
      + ((size_t) 1);
    SM = __alloc_mem(SMs, sizeof(*SM));
    SN = __alloc_mem(SNs, sizeof(*SN));
    __m_memset(SN, 0, SNs, sizeof(*SN));
    __m_memcpy(SM, M, pr, sizeof(*SM));
    __m_memcpy(SN, N, qr, sizeof(*SN));
    shift_left(SN, SNs, (size_t) (-sigma));
  }

  /* We need to return

     (-1)^s * 2^G * (SM + (-1)^(t - s) * SN)

     Let's call

     SM + (-1)^(t - s) * SN

     K

  */
  Ks = SMs;
  if (SNs > Ks) Ks = SNs;
  Ks++;
  K = __alloc_mem(Ks, sizeof(*K));
  __m_memset(K, 0, Ks, sizeof(*K));
  if ((!!s - !!t) != 0) {
    /* K = SM - SN

       We need to be a little careful.

       SN might be greater than SM.

       In that case we need to flip
       SM and SN and flip the sign s.

    */
    SMM = __alloc_mem(Ks, sizeof(*SMM));
    __m_memset(SMM, 0, Ks, sizeof(*SMM));
    __m_memcpy(SMM, SM, SMs, sizeof(*SMM));
    SNN = __alloc_mem(Ks, sizeof(*SNN));
    __m_memset(SNN, 0, Ks, sizeof(*SNN));
    __m_memcpy(SNN, SN, SNs, sizeof(*SNN));
    if (comparison(SMM, SNN, Ks) >= 0) {
      subtraction(K, SMM, Ks, SNN, Ks);
    } else {
      subtraction(K, SNN, Ks, SMM, Ks);
      s = !s;
    }
    __free_mem(SMM);
    __free_mem(SNN);
  } else {
    /* K = SM + SN */
    addition(K, SM, SMs, SN, SNs);
  }
  /* Here, we return

     (-1)^s * 2^G * K

  */
  widefloat_set_from_scaled_integer(z,
                                    s,
                                    G,
                                    K,
                                    Ks);
  __free_mem(SM);
  __free_mem(SN);
  __free_mem(K);
}

/* Return z = -x

   If x is NaN, return NaN
   If x is +Inf, return -Inf
   If x is -Inf, return +Inf
   Otherwise, x is a number,
   return that number with the
   sign flipped.

   CAUTION:

   The precision of z can be different
   from the one of x! Invention of zeros
   or rounding may be required!

*/
void widefloat_neg(widefloat_t *z,
                   const widefloat_t *x) {
  uint64_t *M;

  if (x->fpclass != FPCLASS_NUMBER) {
    switch (x->fpclass) {
    case FPCLASS_NAN:
      z->fpclass = FPCLASS_NAN;
      break;
    case FPCLASS_POS_INF:
      z->fpclass = FPCLASS_NEG_INF;
      break;
    case FPCLASS_NEG_INF:
      z->fpclass = FPCLASS_POS_INF;
      break;
    default:
      z->fpclass = FPCLASS_NAN;
      break;
    }
    z->sign = !x->sign;
    z->exponent = 0;
    __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
    return;
  }

  M = __alloc_mem(x->mantissa_size, sizeof(*M));
  __m_memcpy(M, x->mantissa, x->mantissa_size, sizeof(*M));
  widefloat_set_from_scaled_integer(z,
                                    !x->sign,
                                    (((int64_t) x->exponent) +
                                     ((int64_t) 1) -
                                     (((int64_t) x->mantissa_size) *
                                      ((int64_t) 64)) +
                                     ((int64_t) WIDEFLOAT_OVERHEAD)),
                                    M,
                                    x->mantissa_size);
  __free_mem(M);
}

/* Sets z to (a rounded version of) x */
void widefloat_set(widefloat_t *z,
                   const widefloat_t *x) {
  uint64_t *M;

  if (x->fpclass != FPCLASS_NUMBER) {
    switch (x->fpclass) {
    case FPCLASS_NAN:
      z->fpclass = FPCLASS_NAN;
      break;
    case FPCLASS_POS_INF:
      z->fpclass = FPCLASS_POS_INF;
      break;
    case FPCLASS_NEG_INF:
      z->fpclass = FPCLASS_NEG_INF;
      break;
    default:
      z->fpclass = FPCLASS_NAN;
      break;
    }
    z->sign = x->sign;
    z->exponent = 0;
    __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
    return;
  }

  M = __alloc_mem(x->mantissa_size, sizeof(*M));
  __m_memcpy(M, x->mantissa, x->mantissa_size, sizeof(*M));
  widefloat_set_from_scaled_integer(z,
                                    x->sign,
                                    (((int64_t) x->exponent) +
                                     ((int64_t) 1) -
                                     (((int64_t) x->mantissa_size) *
                                      ((int64_t) 64)) +
                                     ((int64_t) WIDEFLOAT_OVERHEAD)),
                                    M,
                                    x->mantissa_size);
  __free_mem(M);
}

void widefloat_sub(widefloat_t *z,
                   const widefloat_t *x,
                   const widefloat_t *y) {
  widefloat_t ny;

  widefloat_init(&ny, y->mantissa_size);
  widefloat_neg(&ny, y);
  widefloat_add(z, x, &ny);
  widefloat_clear(&ny);
}

/* Compare 2^E * m with 2^F * n where 1 <= m < 2 and 1 <= n < 2. */
static inline int __widefloat_cmp_helper(int32_t E, const uint64_t *m, size_t ms,
                                         int32_t F, const uint64_t *n, size_t ns) {
  int res;
  uint64_t *t;

  if (E < F) return -1;
  if (E > F) return 1;
  if (ms == ns) return comparison(m, n, ms);
  if (ms > ns) {
    t = __alloc_mem(ms, sizeof(*t));
    __m_memset(t, 0, ms, sizeof(*t));
    __m_memcpy(t, n, ns, sizeof(*t));
    shift_left(t, ms, (size_t) ((ms - ns) << 6));
    res = comparison(m, t, ms);
    __free_mem(t);
    return res;
  }
  t = __alloc_mem(ns, sizeof(*t));
  __m_memset(t, 0, ns, sizeof(*t));
  __m_memcpy(t, m, ms, sizeof(*t));
  shift_left(t, ns, (size_t) ((ns - ms) << 6));
  res = comparison(t, n, ns);
  __free_mem(t);
  return res;
}

/* Compares x and y

   If x < y  returns -1
   If x == y return   0
   If x > y  returns  1

   If x or y is NaN, returns 0.

*/
int widefloat_cmp(const widefloat_t *x,
                  const widefloat_t *y) {

  switch (x->fpclass) {
  case FPCLASS_NAN:
    return 0;
    break;
  case FPCLASS_NEG_INF:
    switch (y->fpclass) {
    case FPCLASS_NAN:
      return 0;
      break;
    case FPCLASS_NEG_INF:
      return 0;
      break;
    case FPCLASS_POS_INF:
      return -1;
      break;
    case FPCLASS_NUMBER:
      return -1;
      break;
    default:
      /* Unreachable */
      break;
    }
    /* Unreachable */
    break;
  case FPCLASS_POS_INF:
    switch (y->fpclass) {
    case FPCLASS_NAN:
      return 0;
      break;
    case FPCLASS_NEG_INF:
      return 1;
      break;
    case FPCLASS_POS_INF:
      return 0;
      break;
    case FPCLASS_NUMBER:
      return 1;
      break;
    default:
      /* Unreachable */
      break;
    }
    /* Unreachable */
    break;
  case FPCLASS_NUMBER:
    switch (y->fpclass) {
    case FPCLASS_NAN:
      return 0;
      break;
    case FPCLASS_NEG_INF:
      return 1;
      break;
    case FPCLASS_POS_INF:
      return -1;
      break;
    case FPCLASS_NUMBER:
      if (is_zero(x->mantissa, x->mantissa_size)) {
        if (is_zero(y->mantissa, y->mantissa_size)) {
          return 0;
        } else {
          if (y->sign) {
            /* x = 0, y < 0 */
            return 1;
          } else {
            /* x = 0, y > 0 */
            return -1;
          }
        }
      } else {
        if (is_zero(y->mantissa, y->mantissa_size)) {
          if (x->sign) {
            /* x < 0, y = 0 */
            return -1;
          } else {
            /* x > 0, y = 0 */
            return 1;
          }
        } else {
          /* x != 0, y != 0 */
          if (x->sign) {
            if (y->sign) {
              /* x < 0, y < 0 */
              return __widefloat_cmp_helper(y->exponent, y->mantissa, y->mantissa_size,
                                            x->exponent, x->mantissa, x->mantissa_size);
            } else {
              /* x < 0, y > 0 */
              return -1;
            }
          } else {
            if (y->sign) {
              /* x > 0, y < 0 */
              return 1;
            } else {
              /* x > 0, y > 0 */
              return __widefloat_cmp_helper(x->exponent, x->mantissa, x->mantissa_size,
                                            y->exponent, y->mantissa, y->mantissa_size);
            }
          }
        }
      }
      break;
    default:
      /* Unreachable */
      break;
    }
    /* Unreachable */
    break;
  default:
    /* Unreachable */
    return 0;
    break;
  }
  /* Unreachable */
  return 0;
}

/* If x is NaN or +/-Inf, does nothing and returns -1.
   If x is zero, sets E to 0 and returns 0.
   Otherwise (x is a non-zero number) sets E to

   E = 2 * floor(x->exponent / 2)

   and

   modifies x such that

   2^E * x' = x

   and returns 0

*/
static int __widefloat_rsqrt_helper(int64_t *E, widefloat_t *x) {
  int64_t ee;

  if (x->fpclass != FPCLASS_NUMBER) return -1;
  if (is_zero(x->mantissa, x->mantissa_size)) {
    *E = (int64_t) 0;
    return 0;
  }
  /* Here, x is a non-zero number */
  ee = (int64_t) x->exponent;
  ee >>= 1;
  ee <<= 1;
  x->exponent -= (int32_t) ee;
  *E = ee;
  return 0;
}

/*  Returns 1/sqrt(x)

    If x is NaN, returns NaN.

    If x is +/-Inf, returns 0.

    If x is +/-0, returns +/- Inf.

    If x < 0, return NaN.

    Otherwise returns an approximation to

    1/sqrt(x).

*/
void widefloat_rsqrt(widefloat_t *z, const widefloat_t *x) {
  widefloat_t t;
  int64_t E;
  widefloat_t y;
  widefloat_t three;
  widefloat_t half;
  widefloat_t t1;
  widefloat_t yr, yrp;
  uint64_t m;
  int done;

  if (x->fpclass != FPCLASS_NUMBER) {
    switch (x->fpclass) {
    case FPCLASS_NAN:
      z->sign = 0;
      z->fpclass = FPCLASS_NAN;
      break;
    case FPCLASS_POS_INF:
      z->sign = 0;
      z->fpclass = FPCLASS_NUMBER;
      break;
    case FPCLASS_NEG_INF:
      z->sign = 1;
      z->fpclass = FPCLASS_NUMBER;
      break;
    default:
      z->fpclass = FPCLASS_NAN;
      break;
    }
    z->exponent = 0;
    __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
    return;
  }

  if (is_zero(x->mantissa, x->mantissa_size)) {
    if (x->sign) {
      z->fpclass = FPCLASS_NEG_INF;
    } else {
      z->fpclass = FPCLASS_POS_INF;
    }
    z->sign = x->sign;
    z->exponent = 0;
    __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
    return;
  }

  if (x->sign) {
    z->fpclass = FPCLASS_NAN;
    z->sign = 0;
    z->exponent = 0;
    __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
    return;
  }

  /* Here, x is always a non-zero positive number

     We need to copy x.

  */
  widefloat_init(&t, x->mantissa_size);
  t.fpclass = x->fpclass;
  t.sign = x->sign;
  t.exponent = x->exponent;
  __m_memcpy(t.mantissa, x->mantissa, t.mantissa_size, sizeof(*t.mantissa));
  if (__widefloat_rsqrt_helper(&E, &t) < 0) {
    widefloat_clear(&t);
    return;
  }

  /* Here 1 <= t < 4.

     2^E * t = x

  */
  widefloat_init(&y, z->mantissa_size << 1);
  widefloat_init(&t1, z->mantissa_size << 1);
  widefloat_init(&three, (size_t) 1);
  widefloat_init(&half, (size_t) 1);
  widefloat_init(&yr, z->mantissa_size);
  widefloat_init(&yrp, z->mantissa_size);
  m = (uint64_t) 3;
  widefloat_set_from_scaled_integer(&y, 0, (int64_t) -2, &m, (size_t) 1);
  widefloat_set_from_scaled_integer(&three, 0, (int64_t) 0, &m, (size_t) 1);
  m = (uint64_t) 1;
  widefloat_set_from_scaled_integer(&half, 0, (int64_t) -1, &m, (size_t) 1);
  widefloat_set(&yr, &y);
  done = 0;
  do {
    widefloat_mul(&t1, &y, &y);
    widefloat_mul(&t1, &t1, &t);
    widefloat_sub(&t1, &three, &t1);
    widefloat_mul(&t1, &y, &t1);
    widefloat_mul(&y, &half, &t1);
    widefloat_set(&yrp, &y);
    if (widefloat_cmp(&yr, &yrp) == 0) {
      done = 1;
    } else {
      widefloat_set(&yr, &yrp);
    }
  } while (!done);

  /* We have the result for 1/sqrt(mantissa of x) in yr.

     yr has the same precision as z.

     We need to set z, while integrating the exponent.

  */
  z->fpclass = yr.fpclass;
  z->sign = yr.sign;
  z->exponent = yr.exponent - ((int32_t) (E / ((int64_t) 2)));
  __m_memcpy(z->mantissa, yr.mantissa, z->mantissa_size, sizeof(*(z->mantissa)));

  /* Clear variables */
  widefloat_clear(&t1);
  widefloat_clear(&three);
  widefloat_clear(&half);
  widefloat_clear(&y);
  widefloat_clear(&t);
  widefloat_clear(&yr);
  widefloat_clear(&yrp);
}

/* z = 1/x

   z = 1/sqrt(x) * 1/sqrt(x)

   Special cases:

   For x = NaN  return NaN
   For x = +Inf return +0
   For x = -Inf return -0
   For x = +0   return +Inf
   For x = -0   return -Inf

   Caution: 1/sqrt(x) cannot be computed if
   x is negative. We need to flip the sign twice
   in this case.

*/
void widefloat_rcpr(widefloat_t *z,
                    const widefloat_t *x) {
  widefloat_t t1;
  widefloat_t t2;

  if (x->fpclass != FPCLASS_NUMBER) {
    switch (x->fpclass) {
    case FPCLASS_NAN:
      z->fpclass = FPCLASS_NAN;
      z->sign = 0;
      break;
    case FPCLASS_POS_INF:
      z->fpclass = FPCLASS_NUMBER;
      z->sign = 0;
      break;
    case FPCLASS_NEG_INF:
      z->fpclass = FPCLASS_NUMBER;
      z->sign = 1;
      break;
    default:
      /* Unreachable */
      z->fpclass = FPCLASS_NAN;
      z->sign = 0;
      break;
    }
    z->exponent = 0;
    __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
    return;
  }

  /* Here x is a number.

     Check if x is zero.

  */
  if (is_zero(x->mantissa, x->mantissa_size)) {
    /* x is +/- 0, return +/- Inf */
    if (x->sign) {
      z->fpclass = FPCLASS_NEG_INF;
    } else {
      z->fpclass = FPCLASS_POS_INF;
    }
    z->sign = x->sign;
    z->exponent = 0;
    __m_memset(z->mantissa, 0, z->mantissa_size, sizeof(*(z->mantissa)));
    return;
  }

  /* Here x is a non-zero number

     We copy x into t1. We take the sign off t1.

     We compute t2 = 1/rsqrt(t1). We multiply t2 by itself.

     We put the sign back.

  */
  widefloat_init(&t1, x->mantissa_size);
  widefloat_set(&t1, x);
  t1.sign = 0;
  widefloat_init(&t2, z->mantissa_size + ((size_t) 1));
  widefloat_rsqrt(&t2, &t1);
  widefloat_mul(z, &t2, &t2);
  z->sign = x->sign;
  widefloat_clear(&t1);
  widefloat_clear(&t2);
}


void widefloat_div(widefloat_t *z,
                   const widefloat_t *x,
                   const widefloat_t *y) {
  widefloat_t r;

  widefloat_init(&r, x->mantissa_size + ((size_t) 1));
  widefloat_rcpr(&r, y);
  widefloat_mul(z, x, &r);
  widefloat_clear(&r);
}

/* m = floor(x)

   If x is NaN, does nothing with m and returns -1.
   If x is +/-Inf, does nothing with m and returns -1.
   If x is negative, does nothing with m and returns -1.
   If x is greater than 2^(64 * n) - 1, does nothing with m and returns -1.

   Otherwise, sets m to floor(x) and returns 0.

*/
int widefloat_get_integer(uint64_t *m,
                          size_t n,
                          const widefloat_t *x) {
  int64_t E;
  uint64_t *M;
  uint64_t lzc;
  uint64_t bits_in_M;
  size_t w;
  size_t s;
  int64_t p;

  if (x->fpclass != FPCLASS_NUMBER) return -1;
  if (is_zero(x->mantissa, x->mantissa_size)) {
    __m_memset(m, 0, n, sizeof(*m));
    return 0;
  }
  if (x->sign) return -1;
  /* Here, we need to compute

     floor(2^exponent * 2^(-p + 1) * mantissa)

     with

     p = 64 * mantissa_size - overhead.

     Let us compute

     E = exponent - p + 1

     first.

  */
  E = ((int64_t) x->exponent)
    - ((((int64_t) x->mantissa_size) * ((int64_t) 64)) -
       ((int64_t) WIDEFLOAT_OVERHEAD))
    + ((int64_t) 1);

  /* Here, we need to compute

     floor(2^E * mantissa)

     Two cases:

     If E >= 0: 2^E * mantissa is an integer. The floor goes away.
     If E <  0: 2^E * mantissa is not an integer; we can implement
                floor() by shifting to the right.

  */
  if (E >= ((int64_t) 0)) {
    /* We need to compute 2^E * mantissa by shifting the mantissa to the left.
       We need to be careful as this may overflow.

       We know that

       2^(p - 1) <= mantissa < 2^p

       with

       p = 64 * mantissa_size - overhead.

       We know that we can represent integers in m less than

       2^(64 * n).

       If

       2^(E + p) >= 2^(64 * n), we cannot fit the result.

       E + p >= 64 * n, we cannot fit the result.

    */
    p = (((int64_t) x->mantissa_size) * ((int64_t) 64))
      - ((int64_t) WIDEFLOAT_OVERHEAD);
    if (E + p >= (((int64_t) n) * ((int64_t) 64))) return -1;
    /* Here, the result does fit */
    s = n;
    if (x->mantissa_size > s) s = x->mantissa_size;
    M = __alloc_mem(s, sizeof(*M));
    __m_memset(M, 0, s, sizeof(*M));
    __m_memcpy(M, x->mantissa, x->mantissa_size, sizeof(*M));
    shift_left(M, s, ((size_t) E));
    __m_memcpy(m, M, n, sizeof(*m));
    __free_mem(M);
    return 0;
  }

  /* We need to compute floor(2^E * mantissa) by shifting the mantissa to
     the right by -E.
  */
  M = __alloc_mem(x->mantissa_size, sizeof(*M));
  __m_memcpy(M, x->mantissa, x->mantissa_size, sizeof(*M));
  shift_right(M, x->mantissa_size, (size_t) (-E));

  /* We need to check if M fits into m of size n */
  if (n >= x->mantissa_size) {
    __m_memset(m, 0, n, sizeof(*m));
    __m_memcpy(m, M, n, sizeof(*m));
    __free_mem(M);
    return 0;
  }
  if (is_zero(M, x->mantissa_size)) {
    __m_memset(m, 0, n, sizeof(*m));
    __free_mem(M);
    return 0;
  }
  lzc = leading_zeros(M, x->mantissa_size);
  bits_in_M = (((uint64_t) x->mantissa_size) * ((uint64_t) 64)) - lzc;
  w = (size_t) ((bits_in_M + ((uint64_t) 63)) / ((uint64_t) 64));
  if (w > n) {
    /* Cannot fit the number into m */
    __free_mem(M);
    return -1;
  }
  /* The number does fit */
  __m_memcpy(m, M, n, sizeof(*m));
  __free_mem(M);
  return 0;
}

void widefloat_set_from_safeinteger(widefloat_t *y,
                                    const safeinteger_t *x) {
  uint64_t t;

  if (x == NULL) {
    y->fpclass = FPCLASS_NAN;
    y->sign = 0;
    return;
  }
  if ((x->size == ((size_t) 0)) ||
      (x->data == NULL)) {
    t = (uint64_t) 0;
    widefloat_set_from_integer(y, &t, (size_t) 1);
    return;
  }
  widefloat_set_from_integer(y, x->data, x->size);
  y->sign = !!x->sign;
}

int widefloat_get_safeinteger(safeinteger_t *y,
                              const widefloat_t *x) {
  widefloat_t abs_x;
  size_t s;
  uint64_t E, EE;
  uint64_t *u;

  if (x->fpclass != FPCLASS_NUMBER) return -1;
  if (is_zero(x->mantissa, x->mantissa_size)) {
    safeinteger_set_from_ui(y, (unsigned long long int) 0);
    return 0;
  }
  widefloat_init(&abs_x, x->mantissa_size);
  if (x->sign) {
    widefloat_neg(&abs_x, x);
  } else {
    widefloat_set(&abs_x, x);
  }
  if (abs_x.exponent <= ((int32_t) 0)) {
    s = (size_t) 1;
  } else {
    E = ((uint64_t) abs_x.exponent) + ((uint64_t) 2);
    E >>= 6;
    E++;
    s = (size_t) E;
    EE = (uint64_t) s;
    if (E != EE) {
      widefloat_clear(&abs_x);
      return -1;
    }
    if (s < ((size_t) 1)) s = (size_t) 1;
  }
  u = __alloc_mem(s, sizeof(*u));
  if (widefloat_get_integer(u, s, &abs_x) < 0) {
    __free_mem(u);
    widefloat_clear(&abs_x);
    return -1;
  }
  safeinteger_set_from_integer(y, u, s);
  y->sign = !!x->sign;
  __free_mem(u);
  widefloat_clear(&abs_x);
  return 0;
}
