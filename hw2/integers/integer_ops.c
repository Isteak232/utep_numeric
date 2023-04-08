
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

/* cout * 2^64 + s = a + b */
static inline void __halfadder(uint64_t *cout, uint64_t *s,
                               uint64_t a, uint64_t b) {
  uint64_t ccout, ss;

  ss = a + b;
  ccout = (uint64_t) (ss < a);
  *cout = ccout;
  *s = ss;
}

/* cout * 2^64 + s = a + b + cin */
static inline void __fulladder(uint64_t *cout, uint64_t *s,
                               uint64_t a, uint64_t b, uint64_t cin) {
  uint64_t cc1, cc2, ss1, ss2;

  __halfadder(&cc1, &ss1, a, b);
  __halfadder(&cc2, &ss2, ss1, cin);
  *cout = cc1 + cc2;
  *s = ss2;
}

/* s = (a + b) mod 2^(64 * k)

   a has size m

   b has size n

   k = max(m, n)

*/
void addition(uint64_t *s,
              const uint64_t *a, size_t m,
              const uint64_t *b, size_t n) {
  size_t i;
  uint64_t cin, cout;

  if (m <= n) {
    /* a is shorter than b */
    cin = (uint64_t) 0;
    for (i=0;i<m;i++) {
      __fulladder(&cout, &s[i], a[i], b[i], cin);
      cin = cout;
    }
    /* Invent zeros for a */
    for (;i<n;i++) {
      __halfadder(&cout, &s[i], b[i], cin);
      cin = cout;
    }
  } else {
    /* b is shorter than a */
    cin = (uint64_t) 0;
    for (i=0;i<n;i++) {
      __fulladder(&cout, &s[i], a[i], b[i], cin);
      cin = cout;
    }
    /* Invent zeros for b */
    for (;i<m;i++) {
      __halfadder(&cout, &s[i], a[i], cin);
      cin = cout;
    }
  }
}

/* s = (a - b) mod 2^(64 * k)

   a has size m

   b has size n

   k = max(m, n)

   Subtraction can be implemented as addition.

   a - b = a + (-b) = a + (~b + 1) = a + ~b + 1

   where ~b is a bit flip for each bit in b.

   The extra 1 can be put into the starting carry.

*/
void subtraction(uint64_t *s,
                 const uint64_t *a, size_t m,
                 const uint64_t *b, size_t n) {
  size_t i;
  uint64_t cin, cout;

  if (m <= n) {
    /* a is shorter than b */
    cin = (uint64_t) 1;
    for (i=0;i<m;i++) {
      __fulladder(&cout, &s[i], a[i], ~b[i], cin);
      cin = cout;
    }
    /* Invent zeros for a */
    for (;i<n;i++) {
      __halfadder(&cout, &s[i], ~b[i], cin);
      cin = cout;
    }
  } else {
    /* b is shorter than a */
    cin = (uint64_t) 1;
    for (i=0;i<n;i++) {
      __fulladder(&cout, &s[i], a[i], ~b[i], cin);
      cin = cout;
    }
    /* Invent zeros for b */
    for (;i<m;i++) {
      __fulladder(&cout, &s[i], a[i], ~((uint64_t) 0), cin);
      cin = cout;
    }
  }
}

/* a gets modified such that

   a' = (a * 2^k) mod 2^(64 * n)

*/
void shift_left(uint64_t *a, size_t n, size_t k) {
  size_t w, u, i;
  uint64_t cin, cout;

  /* Shift by 0 bits => nothing to do */
  if (k == ((size_t) 0)) return;

  /* Shift by at least n * 64 bits => set
     result to zero.
  */
  if (k >= (n * ((size_t) 64))) {
    /* Helper function to set memory to zero */
    __m_memset(a, 0, n, sizeof(*a));
    return;
  }

  /* Here we have 1 <= k <= 64 * n - 1

     First "cut" k into w and u such that

     k = w * 64 + u

     Then perform a word shift by w words of
     64 bits.

     Finally perform a bit shift by u bits.

     We have

     k = w * 64 + u

     This means

     u = k - 64 * w

  */
  w = k >> 6;
  u = k - (w << 6);

  /* Do the word shift, if needed */
  if (w >= ((size_t) 1)) {
    for (i=(n-((size_t) 1));i>=w;i--) {
      a[i] = a[i-w];
    }
    for (i=0;i<w;i++) {
      a[i] = (uint64_t) 0;
    }
  }

  /* Do the bit shift, if needed */
  if (u >= ((size_t) 1)) {
    cin = (uint64_t) 0;
    for (i=0;i<n;i++) {
      cout = a[i] >> (64 - u);
      a[i] = (a[i] << u) | cin;
      cin = cout;
    }
  }
}

/* a gets modified such that

   a' = floor(a / 2^k)

*/
void shift_right(uint64_t *a, size_t n, size_t k) {
  size_t w, u, i, l;
  uint64_t cin, cout;

  /* Shift by 0 bits => nothing to do */
  if (k == ((size_t) 0)) return;

  /* Shift by at least n * 64 bits => set
     result to zero.
  */
  if (k >= (n * ((size_t) 64))) {
    /* Helper function to set memory to zero */
    __m_memset(a, 0, n, sizeof(*a));
    return;
  }

  /* Here we have 1 <= k <= 64 * n - 1

     First "cut" k into w and u such that

     k = w * 64 + u

     Then perform a word shift by w words of
     64 bits.

     Finally perform a bit shift by u bits.

     We have

     k = w * 64 + u

     This means

     u = k - 64 * w

  */
  w = k >> 6;
  u = k - (w << 6);

  /* Do the word shift, if needed */
  if (w >= ((size_t) 1)) {
    for (i=0;i<(n-w);i++) {
      a[i] = a[i+w];
    }
    for (i=(n-w);i<n;i++) {
      a[i] = (uint64_t) 0;
    }
  }

  /* Do the bit shift, if needed */
  if (u >= ((size_t) 1)) {
    cin = (uint64_t) 0;
    for (l=n,i=(n-((uint64_t) 1));l>=1;l--,i--) {
      cout = a[i] << (64 - u);
      a[i] = (a[i] >> u) | cin;
      cin = cout;
    }
  }
}

/* a becomes what is in the decimal string mod 2^(64 * n)

   Returns 0 if success
   Returns -1 if failure  (someone tries to convert "cheese")

*/
int convert_from_decimal_string(uint64_t *a, size_t n,
                                const char *str) {
  const char *curr;
  uint64_t digit;
  uint64_t *a_eight;
  uint64_t *a_two;
  uint64_t *t;

  /* Do nothing for the empty string */
  if (str[0] == '\0') return -1;

  /* Allocate memory */
  a_eight = __alloc_mem(n, sizeof(*a_eight));
  a_two = __alloc_mem(n, sizeof(*a_two));
  t = __alloc_mem(n, sizeof(*t));

  /* Set a to zero */
  __m_memset(a, 0, n, sizeof(*a));

  /* Loop over the string */
  for (curr=str; *curr!='\0'; curr++) {
    if (!(('0' <= *curr) &&
          (*curr <= '9'))) {
      /* Free memory */
      __free_mem(a_eight);
      __free_mem(a_two);
      __free_mem(t);

      /* Indicate failure */
      return -1;
    }
    digit = (uint64_t) (((int) *curr) - ((int) '0'));

    /* Multiply a by 10 and add in the digit

       a * 10 = a * 8 + a * 2

    */
    __m_memcpy(a_eight, a, n, sizeof(*a));
    __m_memcpy(a_two, a, n, sizeof(*a));
    shift_left(a_eight, n, 3);
    shift_left(a_two, n, 1);
    addition(t, a_eight, n, a_two, n);
    addition(a, t, n, &digit, 1);
  }

  /* Free memory */
  __free_mem(a_eight);
  __free_mem(a_two);
  __free_mem(t);

  /* Indicate success */
  return 0;
}

/* Returns  -1  if a < b
             0  if a = b
             1  if a > b
*/
int comparison(const uint64_t *a, const uint64_t *b, size_t n) {
  size_t i, k;

  for (i=(n-((size_t) 1)),k=n;k>((size_t) 0);i--,k--) {
    if (a[i] < b[i]) return -1;
    if (a[i] > b[i]) return 1;
  }
  return 0;
}

/* hi * 2^64 + lo = a * b */
static inline void __multiply_digits(uint64_t *hi, uint64_t *lo,
                                     uint64_t a, uint64_t b) {
  uint32_t ah, al, bh, bl;
  uint64_t tah, tal, tbh, tbl;
  uint64_t hh, hl, lh, ll;
  uint32_t llh, lll, lhh, lhl, hlh, hll;
  uint64_t tllh, tlll, tlhh, tlhl, thlh, thll;
  uint32_t t1h, t1l;
  uint64_t t0, t1, t2, tt1h, tt1l, h, l;

  /* Cut into parts */
  ah = (uint32_t) (a >> 32);
  al = (uint32_t) a;
  bh = (uint32_t) (b >> 32);
  bl = (uint32_t) b;
  tah = (uint64_t) ah;
  tal = (uint64_t) al;
  tbh = (uint64_t) bh;
  tbl = (uint64_t) bl;

  /* Partial products */
  hh = tah * tbh;
  hl = tah * tbl;
  lh = tal * tbh;
  ll = tal * tbl;

  /* Cut ll, hl, lh */
  llh = (uint32_t) (ll >> 32);
  lll = (uint32_t) ll;
  hlh = (uint32_t) (hl >> 32);
  hll = (uint32_t) hl;
  lhh = (uint32_t) (lh >> 32);
  lhl = (uint32_t) lh;
  tllh = (uint64_t) llh;
  tlll = (uint64_t) lll;
  thlh = (uint64_t) hlh;
  thll = (uint64_t) hll;
  tlhh = (uint64_t) lhh;
  tlhl = (uint64_t) lhl;

  /* Add parts of same weight */
  t0 = tlll;
  t1 = tllh + thll + tlhl;
  t2 = hh + thlh + tlhh;

  /* Cut t1 into two parts */
  t1h = (uint32_t) (t1 >> 32);
  t1l = (uint32_t) t1;
  tt1h = (uint64_t) t1h;
  tt1l = (uint64_t) t1l;

  /* Add stuff together, reweighting */
  l = t0 + (tt1l << 32);
  h = t2 + tt1h;

  /* Write results back */
  *hi = h;
  *lo = l;
}

/* hi * 2^64 + lo = a * b + c */
static inline void __multiply_and_add(uint64_t *hi, uint64_t *lo,
                                      uint64_t a, uint64_t b, uint64_t c) {
  uint64_t h, l, th, tl, cout;

  __multiply_digits(&th, &tl, a, b);
  __halfadder(&cout, &l, tl, c);
  h = th + cout;
  *hi = h;
  *lo = l;
}

/* p = a * b

   a is on m digits
   b is on 1 digit

   p is on m+1 digits

*/
static inline void __multiplication_rectangular(uint64_t *p,
                                                const uint64_t *a,
                                                size_t m,
                                                uint64_t b) {
  uint64_t cin, cout;
  size_t i;

  /* m is zero: nothing to do */
  if (m == ((size_t) 0)) return;

  /* Loop */
  cin = (uint64_t) 0;
  for (i=0;i<m;i++) {
    __multiply_and_add(&cout, &p[i], a[i], b, cin);
    cin = cout;
  }
  p[i] = cin;
}

static inline unsigned int __floor_log2_size(size_t x) {
  int k;
  size_t t;

  /* If x is zero, return the least value we can return */
  if (x == ((size_t) 0)) return (unsigned int) 0;

  /* If x is one, return 0. */
  if (x == ((size_t) 1)) return (unsigned int) 0;

  /* Here x is at least 2

     There are all kind of smart algorithms.
     We do this with a kind of suboptimal loop.

  */
  for (k=((((int) sizeof(x)) * 8) - 1);
       k>=0;
       k--) {
    t = ((size_t) 1) << ((unsigned int) k);
    if ((x & t) != ((size_t) 0)) {
      return (unsigned int) k;
    }
  }
  /* Unreachable */
  return (unsigned int) 0;
}

/* r = hi * 2^(64 * s) + mi * 2^(64 * t) + lo

   where

   hi has size m,
   mi has size n,
   lo has size p and
   r has size q.

*/
static inline  void __multiplication_helper_sum(uint64_t *r,
                                                size_t q,
                                                const uint64_t *hi,
                                                size_t m, size_t s,
                                                const uint64_t *mi,
                                                size_t n, size_t t,
                                                const uint64_t *lo,
                                                size_t p) {
  uint64_t cin, cout;
  size_t i, z;

  /* Set r to zero. */
  __m_memset(r, 0, q, sizeof(*r));

  /* Copy lo into r */
  z = p; if (q < z) z = q;
  __m_memcpy(r, lo, z, sizeof(*r));

  /* Add in mi */
  cin = (uint64_t) 0;
  if (n+t<q) {
    z = n;
  } else {
    z = q - t;
  }
  for (i=0;i<z;i++) {
    __fulladder(&cout, &r[i+t], r[i+t], mi[i], cin);
    cin = cout;
  }
  if (i+t<q) {
    r[i+t] = cin;
  }

  /* Add in hi */
  cin = (uint64_t) 0;
  if (m+s<q) {
    z = m;
  } else {
    z = q - s;
  }
  for (i=0;i<z;i++) {
    __fulladder(&cout, &r[i+s], r[i+s], hi[i], cin);
    cin = cout;
  }
}

/* p = a * b

   a is on 2^k digits
   b is on 2^k digits

   p is on 2^(k+1) digits

   with

   0 <= k <= 8 * sizeof(size_t) - 1.

*/
static inline void __multiplication_square_karatsuba(uint64_t *p,
                                                     const uint64_t *a,
                                                     const uint64_t *b,
                                                     unsigned int k) {
  unsigned int kprime;
  size_t n;
  uint64_t *ah;
  uint64_t *al;
  uint64_t *bh;
  uint64_t *bl;
  uint64_t *s1;
  uint64_t *s2;
  uint64_t *t1;
  uint64_t *t2;
  uint64_t *t0;
  uint64_t *t2p;
  uint64_t *t0p;
  uint64_t *w;
  int sign_s1, sign_s2;

  /* Check for the base case */
  if (k == ((size_t) 0)) {
    /* Base case: a and b are on 2^0 = 1 digits */
    __multiply_digits(&p[1], &p[0], a[0], b[0]);
    return;
  }

  /* Here, k >= 1.

     We cut a and b into two halves each.

  */
  kprime = k - ((unsigned int) 1);

  /* Allocate memory for ah, al, bh, bl */
  n = ((size_t) 1) << kprime; /* n = 2^(k - 1) */
  ah = __alloc_mem(n, sizeof(*ah));
  al = __alloc_mem(n, sizeof(*al));
  bh = __alloc_mem(n, sizeof(*bh));
  bl = __alloc_mem(n, sizeof(*bl));

  /* Cut a into ah and al, b into bh and bl */
  __m_memcpy(al, a, n, sizeof(*al));
  __m_memcpy(ah, &a[n], n, sizeof(*ah));
  __m_memcpy(bl, b, n, sizeof(*bl));
  __m_memcpy(bh, &b[n], n, sizeof(*bh));

  /* Allocate memory for s1 and s2 */
  s1 = __alloc_mem(n, sizeof(*s1));
  s2 = __alloc_mem(n, sizeof(*s2));

  /* Compute

     s1 = ah - al    resp. al - ah

     s2 = bh - bl    resp. bl - bh

  */
  if (comparison(ah, al, n) >= 0) {
    sign_s1 = 0;
    subtraction(s1, ah, n, al, n);
  } else {
    sign_s1 = 1;
    subtraction(s1, al, n, ah, n);
  }
  if (comparison(bh, bl, n) >= 0) {
    sign_s2 = 0;
    subtraction(s2, bh, n, bl, n);
  } else {
    sign_s2 = 1;
    subtraction(s2, bl, n, bh, n);
  }

  /* Allocate memory for t0, t1, t2, w */
  t0 = __alloc_mem(n, ((size_t) 2) * sizeof(*t0));
  t1 = __alloc_mem(n + ((size_t) 1), ((size_t) 2) * sizeof(*t1));
  t2 = __alloc_mem(n, ((size_t) 2) * sizeof(*t2));
  w = __alloc_mem(n, ((size_t) 2) * sizeof(*w));
  t0p = __alloc_mem(n + ((size_t) 1), ((size_t) 2) * sizeof(*t0p));
  t2p = __alloc_mem(n + ((size_t) 1), ((size_t) 2) * sizeof(*t2p));

  /* Execute the 3 recursive calls */
  __multiplication_square_karatsuba(t0, al, bl, kprime);
  __multiplication_square_karatsuba(t2, ah, bh, kprime);
  __multiplication_square_karatsuba(w, s1, s2, kprime);

  /* Deduce t1 out of w, t0, t2, sign_s1 and sign_s2 */
  __m_memset(t0p, 0, n + ((size_t) 1), ((size_t) 2) * sizeof(*t0p));
  __m_memset(t2p, 0, n + ((size_t) 1), ((size_t) 2) * sizeof(*t2p));
  __m_memcpy(t0p, t0, n, ((size_t) 2) * sizeof(*t0p)); /* t0p = t0 */
  __m_memcpy(t2p, t2, n, ((size_t) 2) * sizeof(*t2p)); /* t2p = t2 */
  addition(t1, t0p, ((size_t) 2) * n + ((size_t) 2), t2p, ((size_t) 2) * n + ((size_t) 2));
  if (sign_s1 + sign_s2 == 1) {
    addition(t1, t1, ((size_t) 2) * n + ((size_t) 2), w, ((size_t) 2) * n);
  } else {
    subtraction(t1, t1, ((size_t) 2) * n + ((size_t) 2), w, ((size_t) 2) * n);
  }

  /* Sum up t2, t1 and t0, all scaled appropriately */
  __multiplication_helper_sum(p, ((size_t) 1) << (k + ((unsigned int) 1)),
                              t2, ((size_t) 2) * n, ((size_t) 1) << k,
                              t1, ((size_t) 2) * n + ((size_t) 2), n,
                              t0, ((size_t) 2) * n);

  /* Free the temporaries */
  __free_mem(ah);
  __free_mem(al);
  __free_mem(bh);
  __free_mem(bl);
  __free_mem(s1);
  __free_mem(s2);
  __free_mem(t0);
  __free_mem(t0p);
  __free_mem(t1);
  __free_mem(t2);
  __free_mem(t2p);
  __free_mem(w);
}

/* Forward declaration */
static inline void __multiplication_square(uint64_t *p,
                                           const uint64_t *a,
                                           const uint64_t *b,
                                           size_t m);

/* p = a * b

   a is on m digits
   b is on m digits

   p is on 2*m digits

   t > m

   t is a power of 2

   The function works only for t > m.

*/
static inline void __multiplication_square_aux1(uint64_t *p,
                                                const uint64_t *a,
                                                const uint64_t *b,
                                                size_t m,
                                                size_t t) {
  uint64_t *aa;
  uint64_t *bb;
  uint64_t *r;

  /* Check the pre-condition */
  if (!(t > m)) return;

  /* Allocate memory */
  aa = __alloc_mem(t, sizeof(*aa));
  bb = __alloc_mem(t, sizeof(*bb));
  r = __alloc_mem(t, ((size_t) 2) * sizeof(*r));

  /* Extend a and b to size t */
  __m_memcpy(aa, a, m, sizeof(*aa));
  __m_memset(&aa[m], 0, (t - m), sizeof(*aa));
  __m_memcpy(bb, b, m, sizeof(*bb));
  __m_memset(&bb[m], 0, (t - m), sizeof(*bb));

  /* Compute r = aa * bb */
  __multiplication_square(r, aa, bb, t);

  /* Compute the 2 * m last digits from r into p */
  __m_memcpy(p, r, (m + m), sizeof(*p));

  /* Free memory */
  __free_mem(aa);
  __free_mem(bb);
  __free_mem(r);
}

/* p = a * b

   a is on m digits
   b is on m digits

   p is on 2*m digits

   m > t

   t is a power of 2

   The function works only for m > t.

   On a 64 bit system, this function is essentially never needed. So
   its performance can be pretty bad.

*/
static inline void __multiplication_square_aux2(uint64_t *p,
                                                const uint64_t *a,
                                                const uint64_t *b,
                                                size_t m,
                                                size_t t) {
  uint64_t *al;
  uint64_t *ah;
  uint64_t *bl;
  uint64_t *bh;
  uint64_t *hh;
  uint64_t *hl;
  uint64_t *lh;
  uint64_t *ll;
  uint64_t *hle;
  uint64_t *lhe;
  uint64_t *hllh;

  /* Check the pre-condition */
  if (!(m > t)) return;

  /* Allocate memory */
  al = __alloc_mem(t, sizeof(*al));
  ah = __alloc_mem(m - t, sizeof(*ah));
  bl = __alloc_mem(t, sizeof(*bl));
  bh = __alloc_mem(m - t, sizeof(*bh));
  hh = __alloc_mem(m - t, ((size_t) 2) * sizeof(*hh));
  hl = __alloc_mem(m, sizeof(*hl));
  lh = __alloc_mem(m, sizeof(*lh));
  ll = __alloc_mem(t, ((size_t) 2) * sizeof(*ll));
  hle = __alloc_mem(m + ((size_t) 1), sizeof(*hle)); /* +1 may overflow */
  lhe = __alloc_mem(m + ((size_t) 1), sizeof(*lhe)); /* +1 may overflow */
  hllh = __alloc_mem(m + ((size_t) 1), sizeof(*hllh)); /* +1 may overflow */

  /* Cut

     a = ah * 2^(64 * t) + al

     and

     b = bh * 2^(64 * t) + bl

  */
  __m_memcpy(al, a, t, sizeof(*al));
  __m_memcpy(ah, &a[t], (m - t), sizeof(*ah));
  __m_memcpy(bl, b, t, sizeof(*bl));
  __m_memcpy(bh, &b[t], (m - t), sizeof(*bh));

  /* Perform the 4 partial products */
  multiplication(hh, ah, m - t, bh, m - t);
  multiplication(hl, ah, m - t, bl, t);
  multiplication(lh, al, t, bh, m - t);
  multiplication(ll, al, t, bl, t);

  /* Add the middle products together */
  __m_memset(hle, 0, (m + ((size_t) 1)), sizeof(*hle));
  __m_memcpy(hle, hl, m, sizeof(*hle));
  __m_memset(lhe, 0, (m + ((size_t) 1)), sizeof(*lhe));
  __m_memcpy(lhe, lh, m, sizeof(*lhe));
  addition(hllh, hle, m + ((size_t) 1), lhe, m + ((size_t) 1));

  /* Put everything back into the result */
  __multiplication_helper_sum(p,
                              m + m,
                              hh, (m - t) + (m - t), t + t,
                              hllh, m + ((size_t) 1), t,
                              ll, t + t);

  /* Free memory */
  __free_mem(al);
  __free_mem(ah);
  __free_mem(bl);
  __free_mem(bh);
  __free_mem(hh);
  __free_mem(hl);
  __free_mem(lh);
  __free_mem(ll);
  __free_mem(hle);
  __free_mem(lhe);
  __free_mem(hllh);
}

/* p = a * b

   a is on m digits
   b is on m digits

   p is on 2*m digits

*/
static inline void __multiplication_square(uint64_t *p,
                                           const uint64_t *a,
                                           const uint64_t *b,
                                           size_t m) {
  unsigned int k;
  size_t t, tt, T;

  /* If m is zero, do nothing */
  if (m == ((size_t) 0)) return;

  /* If m is one, just call the digit multiplication */
  if (m == ((size_t) 1)) {
    __multiply_digits(&p[1], &p[0], a[0], b[0]);
    return;
  }

  /* Here, m >= 2

     Compute

     k = floor(log2(m))

     We have 1 <= k <= 8 * sizeof(size_t) - 1.

  */
  k = __floor_log2_size(m);

  /* If m is equal to 2^k, we can call Karatsuba directly. */
  t = ((size_t) 1) << k;
  if (m == t) {
    __multiplication_square_karatsuba(p, a, b, k);
    return;
  }

  /* Here 2^k = t < m < 2^(k + 1)

     Compute T, the largest power of 2 such that
     32 * T is representable on a size_t.

  */
  T = ((~((size_t) 0)) >> 6) + ((size_t) 1);

  /* If t is less than T, we can compute 2 * t (and 4 * t) without
     overflowing.
  */
  if (t < T) {
    /* We can compute tt = 2 * t, which is a power of 2 and that is
       greater than m.

       We extend a and b to size tt and call ourselves.

       We do all this in a helper function.

    */
    tt = t << 1;
    __multiplication_square_aux1(p, a, b, m, tt);
    return;
  }

  /* Here m > t >= T and T is the largest power of 2 such
     that 2 * T holds on a size_t.

     We need to cut a and b into smaller parts. We do so in
     a helper function.

  */
  __multiplication_square_aux2(p, a, b, m, T);
}

/* p = a * b

   Works only for:

   +  a on m digits
   +  b on n digits
   +  2 <= m
   +  m < n

*/
static inline void __multiplication_aux(uint64_t *p,
                                        const uint64_t *a,
                                        size_t m,
                                        const uint64_t *b,
                                        size_t n) {
  uint64_t *t;
  uint64_t *r;

  /* Handle preconditions */
  if (!(((size_t) 2) <= m)) return;
  if (!(m < n)) return;

  /* Allocate memory */
  t = __alloc_mem(n, sizeof(*t));
  r = __alloc_mem(n, ((size_t) 2) * sizeof(*r));

  /* Copy a into t */
  __m_memcpy(t, a, m, sizeof(*t));

  /* Set upper part of t to zero */
  __m_memset(&t[m], 0, (n - m), sizeof(*t));

  /* Now t and b have the same size n */
  __multiplication_square(r, t, b, n);

  /* Copy low part of r into p */
  __m_memcpy(p, r, (m + n), sizeof(*p));

  /* Free memory */
  __free_mem(t);
  __free_mem(r);
}

/* p = a * b

   a is on m "digits"
   b is on n "digits"

   p must have m + n "digits"

*/
void multiplication(uint64_t *p,
                    const uint64_t *a, size_t m,
                    const uint64_t *b, size_t n) {

  /* If one of the sizes is zero, we do nothing */
  if (m == ((size_t) 0)) return;
  if (n == ((size_t) 0)) return;

  /* If m > n, we flip the arguments */
  if (m > n) {
    multiplication(p, b, n, a, m);
    return;
  }

  /* Here, m <= n

     If m is 1 or n is 1, call a rectangular
     schoolbook multiplication.

  */
  if (m == ((size_t) 1)) {
    __multiplication_rectangular(p, b, n, a[0]);
    return;
  }
  if (n == ((size_t) 1)) {
    __multiplication_rectangular(p, a, m, b[0]);
    return;
  }

  /* Here, 2 <= m <= n.

     If m = n, we call a square multiplication.

  */
  if (m == n) {
    __multiplication_square(p, a, b, m);
    return;
  }

  /* Here, 2 <= m < n.

     We extend a to the size of n.

     We do all that in a helper function.

  */
  __multiplication_aux(p, a, m, b, n);
}


/* Set r = floor(1/10 * 2^(64 * n) */
static inline void __one_tenth(uint64_t *r, size_t n) {
  size_t i;

  /* If n is zero, do nothing */
  if (n == ((size_t) 0)) return;

  /* Set high-level word

     floor(1/10 * 2^64) = 0x1999999999999999

  */
  r[n - ((size_t) 1)] = (uint64_t) 0x1999999999999999ull;

  /* If n is one, we are done */
  if (n == ((size_t) 1)) return;

  /* Here n is at least 2.

     Set the low level words.

     floor(1/10 * 2^(64 * n)) mod 2^64 = 0x9999999999999999

  */
  for (i=0;i<n-((size_t) 1);i++) {
    r[i] = (uint64_t) 0x9999999999999999ull;
  }
}

/* Set

   q = floor(a / 10)

   and

   r = a - 10 * q

   The size of q and a is n.

   r is guaranteed to be 0 <= r <= 9.

*/
void divide_by_ten(uint64_t *q, unsigned int *r, const uint64_t *a, size_t n) {
  uint64_t *one_tenth;
  uint64_t *nine;
  uint64_t *t;
  uint64_t *rr;
  uint64_t *c;
  uint64_t *d;
  uint64_t *e;
  uint64_t one;
  int okay;

  /* If the size is zero, do nothing */
  if (n == ((size_t) 0)) return;

  /* Set one to 1 */
  one = (uint64_t) 1;

  /* Allocate memory for one_tenth on n+2 digits */
  one_tenth = __alloc_mem(n + ((size_t) 2), sizeof(*one_tenth));

  /* Allocate memory for nine on n digits */
  nine = __alloc_mem(n, sizeof(*nine));

  /* Allocate memory for a temporary on n digits */
  rr = __alloc_mem(n, sizeof(*rr));

  /* Allocate memory for a temporary on n digits */
  c = __alloc_mem(n, sizeof(*c));

  /* Allocate memory for a temporary on n digits */
  d = __alloc_mem(n, sizeof(*d));

  /* Allocate memory for a temporary on n digits */
  e = __alloc_mem(n, sizeof(*e));

  /* Allocate memory for a temporary on 2 * n + 2 digits */
  t = __alloc_mem(n + ((size_t) 1), ((size_t) 2) * sizeof(*t));

  /* Load one_tenth = floor(1/10 * 2^(64 * (n + 2))) */
  __one_tenth(one_tenth, n + ((size_t) 2));

  /* Load nine = 9 */
  __m_memset(nine, 0, n, sizeof(*nine));
  nine[0] = (uint64_t) 9;

  /* Multiply a with floor(1/10 * 2^(64 * (n + 2))) */
  multiplication(t, a, n, one_tenth, n + ((size_t) 2));

  /* Divide the temporary by 2^(64 * (n + 2)) */
  __m_memcpy(q, &t[n + ((size_t) 2)], n, sizeof(*q));

  /* Compute, check and correct the remainder */
  okay = 0;
  do {
    /* Multiply q by 10, yielding c */
    __m_memcpy(d, q, n, sizeof(*c));
    __m_memcpy(e, q, n, sizeof(*d));
    shift_left(d, n, 3);
    shift_left(e, n, 1);
    addition(c, d, n, e, n);

    /* We need to compute

       rr = a - c

       We start by checking if c > a.

    */
    if (comparison(c, a, n) > 0) {
      /* c > a. This means 10 * q > a.

         This means q is too great.

         Subtract 1 from q.

      */
      subtraction(c, q, n, &one, (size_t) 1);
      __m_memcpy(q, c, n, sizeof(*q));
      okay = 0;
    } else {
      /* Here, c <= a. Subtract c from a. */
      subtraction(rr, a, n, c, n);

      /* We know that 0 <= rr.

         We need to check if rr <= 9.

      */
      if (comparison(rr, nine, n) > 0) {
        /* Here rr = a - c = a - 10 * q > 9

           This means that q is too small.

           Add 1 to q.

        */
        addition(c, q, n, &one, (size_t) 1);
        __m_memcpy(q, c, n, sizeof(*q));
        okay = 0;
      } else {
        /* Here, we know that

           0 <= rr <= 9.

           The quotient is hence correct.

        */
        okay = 1;
      }
    }
  } while (!okay);

  /* Here, q is already set. Extract r. */
  *r = (unsigned int) rr[0];

  /* Deallocate temporaries */
  __free_mem(one_tenth);
  __free_mem(nine);
  __free_mem(rr);
  __free_mem(c);
  __free_mem(d);
  __free_mem(e);
  __free_mem(t);
}

/* Returns 1 if a is zero. Returns 0 otherwise */
int is_zero(const uint64_t *a, size_t n) {
  size_t i;

  /* Consider empty values to be zero. */
  if (n == ((size_t) 0)) return 1;

  /* Check all digits */
  for (i=0;i<n;i++) {
    if (a[i] != ((uint64_t) 0)) return 0;
  }

  /* We survived the loop. Everything is zero. */
  return 1;
}

/* str becomes the decimal string corresponding to a.

   str needs to have sufficient length.

*/
void convert_to_decimal_string(char *str, const uint64_t *a, size_t n) {
  uint64_t *q;
  uint64_t *t;
  unsigned int r;
  size_t i, k;
  char c;

  /* If n is zero, set str to the empty string */
  if (n == ((size_t) 0)) {
    str[0] = '\0';
    return;
  }

  /* If a is zero, set str to the string "0" */
  if (is_zero(a, n)) {
    str[0] = '0';
    str[1] = '\0';
    return;
  }

  /* a is not zero.

     Allocate space for two temporaries t and q.

  */
  t = __alloc_mem(n, sizeof(*t));
  q = __alloc_mem(n, sizeof(*q));

  /* Copy a into t */
  __m_memcpy(t, a, n, sizeof(*t));

  /* Loop until t is zero. */
  i = (size_t) 0;
  while (!is_zero(t, n)) {
    /* Divide t by 10, put quotient into q,
       remainder into r.
    */
    divide_by_ten(q, &r, t, n);

    /* Copy q into t */
    __m_memcpy(t, q, n, sizeof(*t));

    /* Convert remainder to a digit */
    str[i] = (char) (((int) r) + ((int) '0'));
    i++;
  }
  /* Set end marker in string */
  str[i] = '\0';

  /* Reverse string of length i */
  for (k=0,i--;k<i;k++,i--) {
    c = str[k];
    str[k] = str[i];
    str[i] = c;
  }

  /* Free the temporaries */
  __free_mem(t);
  __free_mem(q);
}

/* Returns the number of leading zero bits in a 64 bit integer.

   If the integer is zero, 64 is returned.

*/
static inline uint64_t __leading_zeros_uint64(uint64_t a) {
  uint64_t res, t;

  if (a == ((uint64_t) 0)) return (uint64_t) 64;

  res = (uint64_t) 0;
  for (t=a;
       (t & (((uint64_t) 1) << 63)) == ((uint64_t) 0);
       t<<=1) {
    res++;
  }
  return res;
}

/* Returns the number of leading zero bits in the integer a of size n.

   If a is zero, 64 * n is returned.
   If n is zero, 0 is returned.

*/
uint64_t leading_zeros(const uint64_t *a, size_t n) {
  uint64_t res;
  size_t i, k;

  if (n == ((size_t) 0)) return (uint64_t) 0;

  res = (uint64_t) 0;
  for (i=n-((size_t) 1),k=n;k>0;k--,i--) {
    if (a[i] != ((uint64_t) 0)) break;
    res += (uint64_t) 64;
  }
  if (k > ((size_t) 0)) {
    res += __leading_zeros_uint64(a[i]);
  }
  return res;
}

int division(uint64_t *q, uint64_t *r,
             const uint64_t *a, size_t m,
             const uint64_t *b, size_t n) {
  widefloat_t afp, bfp, qfp;
  uint64_t *g;
  uint64_t *gb;
  uint64_t *t;
  uint64_t *tt;
  uint64_t *s;
  uint64_t one;

  if (is_zero(b, n)) return -1;

  /* We need to compute

     q = floor(a / b)

     and

     r = a - q * b.

     We get a first guess of q by using floating-point division.

  */
  widefloat_init(&afp, m);
  widefloat_init(&bfp, n);
  widefloat_init(&qfp, m + ((size_t) 1));
  widefloat_set_from_integer(&afp, a, m);
  widefloat_set_from_integer(&bfp, b, n);
  widefloat_div(&qfp, &afp, &bfp);
  g = __alloc_mem(m, sizeof(*g));
  if (widefloat_get_integer(g, m, &qfp) < 0) {
    /* Should never happen */
    __free_mem(g);
    widefloat_clear(&afp);
    widefloat_clear(&bfp);
    widefloat_clear(&qfp);
    return -1;
  }
  widefloat_clear(&afp);
  widefloat_clear(&bfp);
  widefloat_clear(&qfp);

  /* We have g close to floor(a / b).

     We multiply g by b.

  */
  gb = __alloc_mem(m + n, sizeof(*gb));

  /* Allocate a temporary */
  tt = __alloc_mem(m + n, sizeof(*tt));
  __m_memset(tt, 0, m + n, sizeof(*tt));

  /* t = a */
  t = __alloc_mem(m + n, sizeof(*t));
  __m_memset(t, 0, m + n, sizeof(*t));
  __m_memcpy(t, a, m, sizeof(*t));

  /* s = b */
  s = __alloc_mem(m + n, sizeof(*s));
  __m_memset(s, 0, m + n, sizeof(*s));
  __m_memcpy(s, b, n, sizeof(*s));

  /* Check and correct */
  for (;;) {
    multiplication(gb, g, m, b, n);
    /* If gb > t, the guess is too great. We need to subtract 1 off g. */
    if (comparison(gb, t, m + n) > 0) {
      one = (uint64_t) 1;
      subtraction(g, g, m, &one, 1);
    } else {
      /* Here, gb <= t. We can compute t - gb = a - g * b */
      subtraction(tt, t, m + n, gb, m + n);
      /* If tt >= b, the guess is too small. We need to add 1 to g. */
      if (comparison(tt, s, m + n) >= 0) {
        one = (uint64_t) 1;
        addition(g, g, m, &one, 1);
      } else {
        /* The guess is correct. The remainder is in tt. */
        __m_memcpy(q, g, m, sizeof(*q));
        __m_memcpy(r, tt, n, sizeof(*r));
        break;
      }
    }
  }
  __free_mem(g);
  __free_mem(gb);
  __free_mem(t);
  __free_mem(tt);
  __free_mem(s);
  return 0;
}
