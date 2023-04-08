
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
#include "safeinteger_ops.h"
#include "integer_ops.h"


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

static inline void *__realloc_mem(void *ptr, size_t nmemb, size_t size) {
  size_t s;
  void *new_ptr;

  if (ptr == NULL) return __alloc_mem(nmemb, size);
  if (!__try_size_t_multiply(&s, nmemb, size)) {
    fprintf(stderr, "Cannot allocate memory\n");
    exit(1);
  }
  new_ptr = realloc(ptr, s);
  if (s == ((size_t) 0)) return new_ptr;
  if (new_ptr == NULL) {
    fprintf(stderr, "Cannot allocate memory: %s\n", strerror(errno));
    exit(1);
  }
  return new_ptr;
}

static inline void __free_mem(void *ptr) {
  free(ptr);
}

/* Functions to resize a safeinteger's data field/size */

static inline void __safeinteger_resize(safeinteger_t *x, size_t size) {
  if (x == NULL) return;
  if (x->size == ((size_t) 0)) {
    if (size == ((size_t) 0)) return;
    x->data = __alloc_mem(size, sizeof(*(x->data)));
    x->size = size;
    return;
  }
  if (size == ((size_t) 0)) {
    __free_mem(x->data);
    x->data = NULL;
    x->size = (size_t) 0;
    return;
  }
  if (x->size == size) return;
  x->data = __realloc_mem(x->data, size, sizeof(*(x->data)));
  if (size > x->size) {
    __m_memset(&x->data[x->size], 0, size - x->size, sizeof(*(x->data)));
  }
  x->size = size;
}

static inline void __safeinteger_smart_resize(safeinteger_t *x) {
  size_t s, r;
  uint64_t lzc;

  if (x == NULL) return;
  if (x->size == ((size_t) 0)) return;
  if (is_zero(x->data, x->size)) {
    s = (size_t) 1;
  } else {
    lzc = leading_zeros(x->data, x->size);
    r = (size_t) (lzc / ((uint64_t) 64));
    if (r >= x->size) {
      s = (size_t) 1;
    } else {
      s = x->size - r;
    }
  }
  if (s < ((size_t) 1)) s = (size_t) 1;
  if (s >= x->size) return;
  __safeinteger_resize(x, s);
}

/* Actual, exported functions */

void safeinteger_init(safeinteger_t *x) {
  x->sign = 0;
  x->size = (size_t) 0;
  x->data = NULL;
}

void safeinteger_clear(safeinteger_t *x) {
  if (x->data != NULL) {
    __free_mem(x->data);
  }
  x->sign = 0;
  x->size = (size_t) 0;
  x->data = NULL;
}

void safeinteger_set(safeinteger_t *y,
                     const safeinteger_t *x) {
  __safeinteger_resize(y, x->size);
  if ((x->data != NULL) &&
      (y->data != NULL)) {
    __m_memcpy(y->data, x->data, y->size, sizeof(*(y->data)));
  }
  y->sign = x->sign;
  if ((y->data != NULL) &&
      (y->size != ((size_t) 0))) {
    __safeinteger_smart_resize(y);
  }
}

void safeinteger_set_from_ui(safeinteger_t *y,
                             unsigned long long int x) {
  __safeinteger_resize(y, ((size_t) 1));
  y->sign = 0;
  y->data[0] = (uint64_t) x;
  __safeinteger_smart_resize(y);
}

void safeinteger_set_from_si(safeinteger_t *y,
                             long long int x) {
  unsigned long long int xabs;
  int sign;

  if (x >= ((long long int) 0)) {
    xabs = (unsigned long long int) x;
    sign = 0;
  } else {
    xabs = ((unsigned long long int) (-(x + ((long long int) 1))))
      + ((unsigned long long int) 1);
    sign = 1;
  }
  safeinteger_set_from_ui(y, xabs);
  y->sign = sign;
}

int safeinteger_get_ui(unsigned long long int *y,
                       const safeinteger_t *x) {
  size_t s, r;
  uint64_t lzc;
  uint64_t v, vvv;
  unsigned long long int vv;

  if (x->size == ((size_t) 0)) return -1;
  if (x->data == NULL) return -1;
  if (x->sign) return -1;
  if (is_zero(x->data, x->size)) {
    s = (size_t) 1;
  } else {
    lzc = leading_zeros(x->data, x->size);
    r = (size_t) (lzc / ((uint64_t) 64));
    if (r >= x->size) {
      s = (size_t) 1;
    } else {
      s = x->size - r;
    }
  }
  if (s < ((size_t) 1)) s = (size_t) 1;
  if (s >= x->size) return -1;
  if (s > ((size_t) 1)) return -1;
  v = x->data[0];
  vv = (unsigned long long int) v;
  vvv = (uint64_t) vv;
  if (v != vvv) return -1;
  *y = vv;
  return 0;
}

int safeinteger_get_si(long long int *y,
                       const safeinteger_t *x) {
  size_t s, r;
  uint64_t lzc;
  uint64_t v, vvv;
  unsigned long long int vv, vvvv;
  long long int t;

  if (x->size == ((size_t) 0)) return -1;
  if (x->data == NULL) return -1;
  if (is_zero(x->data, x->size)) {
    s = (size_t) 1;
  } else {
    lzc = leading_zeros(x->data, x->size);
    r = (size_t) (lzc / ((uint64_t) 64));
    if (r >= x->size) {
      s = (size_t) 1;
    } else {
      s = x->size - r;
    }
  }
  if (s < ((size_t) 1)) s = (size_t) 1;
  if (s >= x->size) return -1;
  if (s > ((size_t) 1)) return -1;
  v = x->data[0];
  vv = (unsigned long long int) v;
  vvv = (uint64_t) vv;
  if (v != vvv) return 0;
  if (vv == ((unsigned long long int) 0)) {
    t = (long long int) 0;
  } else {
    vv--;
    t = (long long int) vv;
    if (t < ((long long int) 0)) return -1;
    vvvv = (unsigned long long int) t;
    if (vv != vvvv) return -1;
    if (x->sign) {
      if (t != ((long long int) 0)) {
        t = -t;
        if (t >= ((long long int) 0)) return -1;
        t--;
        if (t >= ((long long int) 0)) return -1;
      }
    } else {
      t++;
      if (t <= ((long long int) 0)) return -1;
    }
  }
  *y = t;
  return 0;
}

void safeinteger_neg(safeinteger_t *y,
                     const safeinteger_t *x) {
  if (x->size == ((size_t) 0)) return;
  if (x->data == NULL) return;
  __safeinteger_resize(y, x->size);
  y->sign = !x->sign;
  __m_memcpy(y->data, x->data, y->size, sizeof(*(y->data)));
  __safeinteger_smart_resize(y);
}

void safeinteger_add(safeinteger_t *z,
                     const safeinteger_t *x,
                     const safeinteger_t *y) {
  size_t s;
  uint64_t *a;
  uint64_t *b;

  if (x->size == ((size_t) 0)) return;
  if (x->data == NULL) return;
  if (y->size == ((size_t) 0)) return;
  if (y->data == NULL) return;

  s = x->size;
  if (y->size > s) s = y->size;
  s++;
  __safeinteger_resize(z, s);
  a = __alloc_mem(s, sizeof(*a));
  b = __alloc_mem(s, sizeof(*b));
  __m_memset(a, 0, s, sizeof(*a));
  __m_memset(b, 0, s, sizeof(*b));
  __m_memset(z->data, 0, z->size, sizeof(*(z->data)));
  __m_memcpy(a, x->data, x->size, sizeof(*a));
  __m_memcpy(b, y->data, y->size, sizeof(*b));
  if ((!!x->sign) == (!!y->sign)) {
    z->sign = x->sign;
    addition(z->data, a, s, b, s);
  } else {
    if (comparison(a, b, s) >= 0) {
      /* 4 + (-3) = 4 - 3  or  -4 + 3 = -(4 - 3) */
      z->sign = x->sign;
      subtraction(z->data, a, s, b, s);
    } else {
      /* 3 + (-4) = -(4 - 3)  or  -3 + 4 = 4 - 3 */
      z->sign = y->sign;
      subtraction(z->data, b, s, a, s);
    }
  }
  __free_mem(a);
  __free_mem(b);
  __safeinteger_smart_resize(z);
}

void safeinteger_sub(safeinteger_t *z,
                     const safeinteger_t *x,
                     const safeinteger_t *y) {
  safeinteger_t t;

  safeinteger_init(&t);
  safeinteger_neg(&t, y);
  safeinteger_add(z, x, &t);
  safeinteger_clear(&t);
}

void safeinteger_mul(safeinteger_t *z,
                     const safeinteger_t *x,
                     const safeinteger_t *y) {
  size_t s;

  if (x->size == ((size_t) 0)) return;
  if (x->data == NULL) return;
  if (y->size == ((size_t) 0)) return;
  if (y->data == NULL) return;
  s = x->size + y->size;
  if (s < x->size) return;
  __safeinteger_resize(z, s);
  multiplication(z->data, x->data, x->size, y->data, y->size);
  z->sign = !!((!!x->sign) ^ (!!y->sign));
  __safeinteger_smart_resize(z);
}

int safeinteger_div(safeinteger_t *quot,
                    safeinteger_t *rema,
                    const safeinteger_t *x,
                    const safeinteger_t *y) {
  size_t qs, rs;

  if (x->size == ((size_t) 0)) return -1;
  if (x->data == NULL) return -1;
  if (y->size == ((size_t) 0)) return -1;
  if (y->data == NULL) return -1;
  if (is_zero(y->data, y->size)) return -1;
  if (is_zero(x->data, x->size)) {
    safeinteger_set(quot, x);
    safeinteger_set(rema, x);
    return 1;
  }
  qs = quot->size;
  rs = rema->size;
  __safeinteger_resize(quot, x->size);
  __safeinteger_resize(rema, y->size);
  if (division(quot->data, rema->data,
                x->data, x->size,
                y->data, y->size) < 0) {
    __safeinteger_resize(quot, qs);
    __safeinteger_resize(rema, rs);
    return -1;
  }
  quot->sign = !!((!!x->sign) ^ (!!y->sign));
  rema->sign = !!x->sign;
  __safeinteger_smart_resize(quot);
  __safeinteger_smart_resize(rema);
  return 0;
}

void safeinteger_shift_left(safeinteger_t *z,
                            const safeinteger_t *x,
                            size_t y) {
  size_t s, e;

  if (y == ((size_t) 0)) {
    safeinteger_set(z, x);
    return;
  }
  if (x->size == ((size_t) 0)) return;
  if (x->data == NULL) return;
  e = (y / ((size_t) 64)) + ((size_t) 1);
  s = x->size + e;
  if (s < x->size) return;
  __safeinteger_resize(z, s);
  z->sign = x->sign;
  __m_memset(z->data, 0, z->size, sizeof(*(z->data)));
  __m_memcpy(z->data, x->data, x->size, sizeof(*(z->data)));
  shift_left(z->data, z->size, y);
  __safeinteger_smart_resize(z);
}

void safeinteger_shift_right(safeinteger_t *z,
                             const safeinteger_t *x,
                             size_t y) {
  uint64_t one;

  if (y == ((size_t) 0)) {
    safeinteger_set(z, x);
    return;
  }
  if (x->size == ((size_t) 0)) return;
  if (x->data == NULL) return;
  __safeinteger_resize(z, x->size);
  z->sign = x->sign;
  if (is_zero(x->data, x->size)) {
    __m_memcpy(z->data, x->data, z->size, sizeof(*(z->data)));
  } else {
    if (x->sign) {
      one = (uint64_t) 1;
      subtraction(z->data, x->data, x->size, &one, ((size_t) 1));
      shift_right(z->data, z->size, y);
      addition(z->data, z->data, z->size, &one, ((size_t) 1));
    } else {
      __m_memcpy(z->data, x->data, z->size, sizeof(*(z->data)));
      shift_right(z->data, z->size, y);
    }
  }
  __safeinteger_smart_resize(z);
}

int safeinteger_cmp(const safeinteger_t *x,
                    const safeinteger_t *y) {
  size_t s;
  uint64_t *a;
  uint64_t *b;
  int res;

  if (x->size == ((size_t) 0)) return 0;
  if (x->data == NULL) return 0;
  if (y->size == ((size_t) 0)) return 0;
  if (y->data == NULL) return 0;
  if ((!!x->sign) != (!!y->sign)) {
    if (is_zero(x->data, x->size) &&
        is_zero(y->data, y->size)) return 0;
    if (x->sign) return -1;
    return 1;
  }
  s = x->size;
  if (y->size > s) s = y->size;
  a = __alloc_mem(s, sizeof(*a));
  b = __alloc_mem(s, sizeof(*b));
  __m_memset(a, 0, s, sizeof(*a));
  __m_memset(b, 0, s, sizeof(*b));
  __m_memcpy(a, x->data, x->size, sizeof(*a));
  __m_memcpy(b, y->data, y->size, sizeof(*b));
  res = comparison(a, b, s);
  if (x->sign) {
    res = -res;
  }
  __free_mem(a);
  __free_mem(b);
  return res;
}

int safeinteger_from_string(safeinteger_t *x,
                            const char *str) {
  size_t len, xs, s;
  const char *abs_str;

  if (str == NULL) return -1;
  len = strlen(str);
  if (len == ((size_t) 0)) return -1;
  s = len << 2u;
  if (len != (s >> 2u)) return -1;
  xs = x->size;
  __safeinteger_resize(x, s);
  abs_str = str;
  if (str[0] == '-') {
    x->sign = 1;
    abs_str++;
    if (len == ((size_t) 1)) {
      __safeinteger_resize(x, xs);
      return -1;
    }
  } else {
    x->sign = 0;
  }
  if (convert_from_decimal_string(x->data, x->size,
                                  abs_str) < 0) {
    __safeinteger_resize(x, xs);
    return -1;
  }
  __safeinteger_smart_resize(x);
  return 0;
}

char *safeinteger_to_string(const safeinteger_t *x) {
  char *str;
  size_t l;
  char *abs_str;

  if (x == NULL) return NULL;
  if ((x->size == ((size_t) 0)) ||
      (x->data == NULL)) {
    str = __alloc_mem(((size_t) 2), sizeof(*str));
    str[0] = '0';
    str[1] = '\0';
    return str;
  }
  if (is_zero(x->data, x->size)) {
    str = __alloc_mem(((size_t) 2), sizeof(*str));
    str[0] = '0';
    str[1] = '\0';
    return str;
  }
  l = (x->size >> 1u) + ((size_t) 2);
  str = __alloc_mem(l, sizeof(*str));
  __m_memset(str, (int) '\0', l, sizeof(*str));
  abs_str = str;
  if (x->sign) {
    str[0] = '-';
    abs_str++;
  }
  convert_to_decimal_string(abs_str, x->data, x->size);
  l = strlen(str);
  l++;
  str = __realloc_mem(str, l, sizeof(*str));
  return str;
}

void safeinteger_to_array(int *sign,
                          uint32_t **abs,
                          size_t *size,
                          const safeinteger_t *x) {
  uint32_t *arr;
  size_t i;
  uint32_t h, l;

  if ((x->size == ((size_t) 0)) ||
      (x->data == NULL)) {
    *sign = 0;
    *size = (size_t) 1;
    *abs = __alloc_mem(*size, sizeof(**abs));
    **abs = (uint32_t) 0;
    return;
  }
  *sign = !!x->sign;
  *size = x->size << 1u;
  *abs = __alloc_mem(*size, sizeof(**abs));
  arr = *abs;
  for (i=(size_t) 0;i<x->size;i++) {
    l = (uint32_t) x->data[i];
    h = (uint32_t) (x->data[i] >> 32u);
    arr[(i << 1u) + ((size_t) 0)] = l;
    arr[(i << 1u) + ((size_t) 1)] = h;
  }
}

void safeinteger_from_array(safeinteger_t *x,
                            int sign,
                            uint32_t *abs,
                            size_t size) {
  size_t s, i, k;
  uint32_t h, l;
  uint64_t t, hh, ll;

  if (size == ((size_t) 0)) {
    safeinteger_set_from_ui(x, ((unsigned long long int) 0));
    return;
  }
  s = (size >> 1u) + ((size_t) 1);
  __safeinteger_resize(x, s);
  x->sign = !!sign;
  for (i=(size_t) 0, k=(size_t) 0; i<size; i+=(size_t) 2, k++) {
    h = abs[i + ((size_t) 1)];
    l = abs[i + ((size_t) 0)];
    hh = (uint64_t) h;
    ll = (uint64_t) l;
    t = (hh << 32u) | ll;
    x->data[k] = t;
  }
  i -= (size_t) 1;
  if (i<size) {
    l = abs[i + ((size_t) 0)];
    ll = (uint64_t) l;
    x->data[k] = ll;
  }
  __safeinteger_smart_resize(x);
}

void safeinteger_set_from_integer(safeinteger_t *y,
                                  uint64_t *x,
                                  size_t size) {
  if (size == ((size_t) 0)) {
    safeinteger_set_from_ui(y, ((unsigned long long int) 0));
    return;
  }
  __safeinteger_resize(y, size);
  __m_memcpy(y->data, x, y->size, sizeof(*(y->data)));
  y->sign = 0;
  __safeinteger_smart_resize(y);
}
