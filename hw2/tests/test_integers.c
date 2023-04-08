 
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "utepnum.h"
#include <string.h>

static int convert_to_size(size_t *res, const char *str) {
  unsigned long long int t, ttt;
  size_t tt;
  char *end;

  if (str == NULL) return -1;
  if (*str == '\0') return -1;
  t = strtoull(str, &end, 0);
  if (*end != '\0') return -1;
  tt = (size_t) t;
  ttt = (unsigned long long int) tt;
  if (ttt != t) return -1;
  *res = tt;
  return 0;
}

static void print_array(const char *str, const uint64_t *a, size_t n) {
  char *sep = "";
  size_t i, k;

  printf("%s[", str);
  if (n == ((size_t) 0)) {
    printf("]\n");
    return;
  }

  for (k=n,i=n-((size_t) 1);k>0;k--,i--) {
    printf("%s%llu", sep, (unsigned long long int) a[i]);
    sep = ",";
  }
  printf("]\n");
}

int test_integers(size_t m, size_t n, const char *str1, const char *str2) {
  size_t q = (m > n) ? m : n;
  size_t r = m + n;
  uint64_t a[m];
  uint64_t b[n];
  uint64_t sum[q];
  uint64_t diff[q];
  uint64_t prod[r];
  uint64_t quot[m];
  uint64_t rema[n];
  uint64_t sl[m];
  uint64_t sr[m];
  char     strres1[r * 64 + 1];
  char     strres2[r * 64 + 1];
  size_t   sigma;


  /* Convert the two strings str1 and str2 */
  if (convert_from_decimal_string(a, m, str1) < 0) return -1;
  if (convert_from_decimal_string(b, n, str2) < 0) return -1;

  /* Display the two arrays */
  print_array("a = ", a, m);
  print_array("b = ", b, n);

  /* Convert a and b back to decimal and print */
  convert_to_decimal_string(strres1, a, m);
  convert_to_decimal_string(strres2, b, n);

  /* Display a and b reconverted to string */
  printf("a = \"%s\"\n", strres1);
  printf("b = \"%s\"\n", strres2);

  /* Test addition and subtraction */
  addition(sum, a, m, b, n);
  subtraction(diff, a, m, b, n);

  /* Display the arrays for sum and difference */
  print_array("sum  = ", sum, q);
  print_array("diff = ", diff, q);

  /* Convert the sum and the difference to decimal and print */
  convert_to_decimal_string(strres1, sum, q);
  convert_to_decimal_string(strres2, diff, q);
  printf("sum  = \"%s\"\n", strres1);
  printf("diff = \"%s\"\n", strres2);

  /* Test multiplication */
  multiplication(prod, a, m, b, n);

  /* Convert the product to decimal and print */
  convert_to_decimal_string(strres1, prod, r);
  printf("prod = \"%s\"\n", strres1);

  /* Test division */
  if (division(quot, rema, a, m, b, n) == 0) {
    printf("Division has signaled success.\n");
  } else {
    printf("Division has signaled failure.\n");
  }

  /* Convert the quotient and the remainder to decimal and print */
  convert_to_decimal_string(strres1, quot, m);
  convert_to_decimal_string(strres2, rema, n);
  printf("quot = \"%s\"\n", strres1);
  printf("rem  = \"%s\"\n", strres2);

  memcpy(sl, a, sizeof(sl));
  memcpy(sr, a, sizeof(sr));

  sigma = ((size_t) b[0]) % ((size_t) 1000);
  shift_left(sl, m, sigma);
  shift_right(sr, m, sigma); /* TODO: CHECK */

  /* Convert the shifted numbers and print */
  convert_to_decimal_string(strres1, sl, m);
  convert_to_decimal_string(strres2, sr, m);
  printf("shift amount: %zu\n", sigma);
  printf("sl = \"%s\"\n", strres1);
  printf("sr = \"%s\"\n", strres2);


  /* Signal success */
  return 0;
}

int main(int argc, char **argv) {
  size_t m, n;

  /* Check if we have at least 5 arguments */
  if (argc < 5) return 1;

  /* Convert the first two arguments to size_t */
  if (convert_to_size(&m, argv[1]) < 0) return 1;
  if (convert_to_size(&n, argv[2]) < 0) return 1;

  /* Check that none of m or n is zero */
  if (m == ((size_t) 0)) return 1;
  if (n == ((size_t) 0)) return 1;

  /* Run the actual test function */
  if (test_integers(m, n, argv[3], argv[4]) < 0) return 1;

  /* Signal success */
  return 0;
}
