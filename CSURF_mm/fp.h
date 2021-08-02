#ifndef FP_H
#define FP_H

#include <stdbool.h>

#include "params.h"

extern const fp fp_0;
extern const fp fp_1;
extern const fp base_a;
extern const fp base_a_inv;
extern const fp fp_2;
extern const fp fp_3;
extern const fp fp_4;
extern const fp fp_8;
extern const fp fp_4inv;

void fp_set(fp *x, uint64_t y);
void fp_cswap(fp *x, fp *y, bool c);
void fp_copy(fp *c, fp *a);


void fp_enc(fp *x, uint const *y); /* encode to Montgomery representation */
void fp_dec(uint *x, fp const *y); /* decode from Montgomery representation */

void fp_add2(fp *x, fp const *y);
void fp_sub2(fp *x, fp const *y);
void fp_mul2(fp *x, fp const *y);

void fp_add3(fp *x, fp const *y, fp const *z);
void fp_sub3(fp *x, fp const *y, fp const *z);
void fp_mul3(fp *x, fp const *y, fp const *z);


void fp_sq1(fp *x);
void fp_sq2(fp *x, fp const *y);
void fp_inv(fp *x);
void fp_sqrt(fp *x);
void fp_triexp(fp *x);
bool fp_issquare(fp *x); /* destroys input! */

void fp_random(fp *x);

#endif