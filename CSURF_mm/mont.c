
#include <assert.h>

#include "params.h"
#include "steps.h"
#include "uint.h"
#include "fp.h"
#include "mont.h"
#include "poly.h"
#include "swindows.h"
#include <stdio.h>
#include <string.h>

#define BITS 512
void xDBLADD(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A24)
{
    fp tmp0, tmp1, tmp2;        //requires precomputation of A24=(A+2C:4C)

    fp_add3(&tmp0, &P->x, &P->z);
    fp_sub3(&tmp1, &P->x, &P->z);
    fp_sq2(&R->x, &tmp0);
    fp_sub3(&tmp2, &Q->x, &Q->z);
    fp_add3(&S->x, &Q->x, &Q->z);
    fp_mul2(&tmp0, &tmp2);
    fp_sq2(&R->z, &tmp1);
    fp_mul2(&tmp1, &S->x);
    fp_sub3(&tmp2, &R->x, &R->z);
    fp_mul2(&R->z, &A24->z);
    fp_mul2(&R->x, &R->z);
    fp_mul3(&S->x, &A24->x, &tmp2);
    fp_sub3(&S->z, &tmp0, &tmp1);
    fp_add2(&R->z, &S->x);
    fp_add3(&S->x, &tmp0, &tmp1);
    fp_mul2(&R->z, &tmp2);
    fp_sq1(&S->z);
    fp_sq1(&S->x);
    fp_mul2(&S->z, &PQ->x);
    fp_mul2(&S->x, &PQ->z);
}


void xDBLADD_min(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A)
{
    fp a, b, c, d, e, f;

    fp_sq2(&a, &P->x);                  // a = X^2
    fp_sq2(&b, &P->z);                  // b = Z^2
    fp_add3(&d, &a, &b);                // d = X^2+Z^2
    fp_sq2(&c, &d);                     // c = (X^2+Z^2)^2

    fp_add3(&e, &P->x, &P->z);          // e = X+Z
    fp_sq1(&e);                         // e = (X+Z)^2
    fp_sub2(&e, &d);                    // e = 2XZ


    fp_sub3(&a, &a, &b);                  // a = X^2-Z^2
    fp_mul2(&a, &A->z);                   // a = c(X^2-Z^2)
    fp_mul3(&b, &A->x, &e);               // b = 2AXZ

    fp_add3(&a, &a, &a);                  // a=2C(X^2-Z^2)
    fp_add2(&a,&b);                       // a=2C(X^2-Z^2)+2AXZ
    // b, d
    fp_mul3(&b, &P->x, &Q->x);            // b=XX_i
    fp_mul3(&d, &P->z, &Q->z);            // d=ZZ_i
    fp_add2(&b,&d);                       // b=XX_i+ZZ_i
    fp_sq1(&b);                           // b=(XX_i+ZZ_i)^2

    fp_mul3(&d, &P->z, &Q->x);          //d=ZX_i
    fp_mul3(&f, &P->x, &Q->z);          //f=XZ_i
    fp_sub2(&d,&f);
    fp_sq1(&d); 
    fp_mul3(&R->z, &a, &e);             //
    fp_mul3(&R->x, &A->z, &c);            // Qx = c(X^2+Z^2)^2
    fp_mul3(&S->x, &PQ->z, &b);
    fp_mul3(&S->z, &PQ->x, &d);
}



void xDBLADD_min_affine(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A)
{
    fp a, b, c, d, e, f;

    fp_sq2(&a, &P->x);                  // a = X^2
    fp_sq2(&b, &P->z);                  // b = Z^2
    fp_add3(&d, &a, &b);                // d = X^2+Z^2
    fp_sub3(&c, &a, &b);

    fp_mul3(&a, &P->x, &P->z);          // a = XZ
    fp_mul3(&b, &a, &A->x);             // b = AXZ
    fp_add3(&a, &a, &a);
    fp_add3(&a, &a, &a);                // a = 4XZ

    fp_add3(&c, &c, &b);                // c = X^2-Z^2+aXZ


    // b, d
    fp_mul3(&b, &P->x, &Q->x);            // b=XX_i
    fp_mul3(&e, &P->z, &Q->z);            // d=ZZ_i
    fp_add2(&b,&e);                       // b=XX_i+ZZ_i
    fp_sq1(&b);                           // b=(XX_i+ZZ_i)^2

    fp_mul3(&e, &P->z, &Q->x);          //d=ZX_i
    fp_mul3(&f, &P->x, &Q->z);          //f=XZ_i
    fp_sub2(&e,&f);
    fp_sq1(&e); 
    fp_mul3(&R->z, &c, &a);
    fp_mul3(&R->x, &d, &d);             // x = (X^2+Z^2)^2


    fp_mul3(&S->x, &PQ->z, &b);
    fp_mul3(&S->z, &PQ->x, &e);
}



void xDBL(proj *Q, proj const *A, proj const *P)
{
    fp a, b, c;
    fp_add3(&a, &P->x, &P->z);
    fp_sq1(&a);
    fp_sub3(&b, &P->x, &P->z);
    fp_sq1(&b);
    fp_sub3(&c, &a, &b);
    fp_add2(&b, &b); fp_add2(&b, &b); /* multiplication by 4 */
    fp_mul2(&b, &A->z);
    fp_mul3(&Q->x, &a, &b);
    fp_add3(&a, &A->z, &A->z); /* multiplication by 2 */
    fp_add2(&a, &A->x);
    fp_mul2(&a, &c);
    fp_add2(&a, &b);
    fp_mul3(&Q->z, &a, &c);
}


void xDBL_min(proj *Q, proj const *A, proj const *P)
{
    fp a, b, c, d, e;
    fp_sq2(&a, &P->x);                  // a = X^2
    fp_sq2(&b, &P->z);                  // b = Z^2
    fp_add3(&d, &a, &b);                // d = X^2+Z^2
    fp_sq2(&c, &d);                     // c = (X^2+Z^2)^2
    fp_add3(&e, &P->x, &P->z);          // e = X+Z   
    fp_sq1(&e);                         // e = (X+Z)^2
    fp_sub2(&e, &d);                    // e = 2XZ
    fp_sub3(&a, &a, &b);                  // a = X^2-Z^2
    fp_mul2(&a, &A->z);                   // a = c(X^2-Z^2)    
    fp_mul3(&b, &A->x, &e);               // b = 2AXZ
    fp_mul3(&Q->x, &A->z, &c);            // Qx = c(X^2+Z^2)^2
    fp_add3(&a, &a, &a);                  // a = 2c(X^2-Z^2)
    fp_add2(&a,&b);                       // a = 2c(X^2-Z^2)+2AXZ
    fp_mul3(&Q->z, &a, &e);               // Qz = 2XZ(2c(X^2-Z^2)+2AXZ)

}

void xADD(proj *S, proj const *P, proj const *Q, proj const *PQ)
{
    fp a, b, c, d;
    fp_add3(&a, &P->x, &P->z);
    fp_sub3(&b, &P->x, &P->z);
    fp_add3(&c, &Q->x, &Q->z);
    fp_sub3(&d, &Q->x, &Q->z);
    fp_mul2(&a, &d);
    fp_mul2(&b, &c);
    fp_add3(&c, &a, &b);
    fp_sub3(&d, &a, &b);
    fp_sq1(&c);
    fp_sq1(&d);
    fp_mul3(&S->x, &PQ->z, &c);
    fp_mul3(&S->z, &PQ->x, &d);
}


void xADD_min(proj *S, proj const *P, proj const *Q, proj const *PQ)
{
    fp a, b, c, d;
    fp_mul3(&a, &P->x, &Q->x);          //a=XX_i
    fp_mul3(&b, &P->z, &Q->z);          //b=ZZ_i
    fp_mul3(&c, &P->z, &Q->x);          //c=ZX_i
    fp_mul3(&d, &P->x, &Q->z);          //d=XZ_i
    fp_add2(&a,&b);
    fp_sub2(&c,&d);
    fp_sq1(&a);                         //n
    fp_sq1(&c);
    fp_mul3(&S->x, &PQ->z, &a);
    fp_mul3(&S->z, &PQ->x, &c);
}




/* Montgomery ladder. */
/* P must not be the unique point of order 2. */
/* not constant-time! */
void xMUL(proj *Q, proj const *A, proj const *P, uint const *k)
{
    proj R = *P;
    proj A24;
    const proj Pcopy = *P; /* in case Q = P */

    Q->x = fp_1;
    Q->z = fp_0;

    fp_add3(&A24.x, &A->z, &A->z);    //precomputation of A24=(A+2C:4C)
    fp_add3(&A24.z, &A24.x, &A24.x);
    fp_add2(&A24.x, &A->x);

    unsigned long i = BITS;
    while (--i && !uint_bit(k, i));

    do {

        bool bit = uint_bit(k, i);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);

        xDBLADD(Q, &R, Q, &R, &Pcopy, &A24);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);

    } while (i--);
}


void xMUL_min(proj *Q, proj const *A, proj const *P, uint const *k)
{
    proj R = *P;
    const proj Pcopy = *P; /* in case Q = P */

    Q->x = fp_1;
    Q->z = fp_0;

    unsigned long i = BITS;
    while (--i && !uint_bit(k, i));

    do {

        bool bit = uint_bit(k, i);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);

        xDBLADD_min(Q, &R, Q, &R, &Pcopy, A);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);

    } while (i--);
}


void xMUL_min_affine(proj *Q, proj const *A, proj const *P, uint const *k)
{
    proj R = *P;
    const proj Pcopy = *P; /* in case Q = P */

    Q->x = fp_1;
    Q->z = fp_0;

    unsigned long i = BITS;
    while (--i && !uint_bit(k, i));

    do {

        bool bit = uint_bit(k, i);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);

        xDBLADD_min_affine(Q, &R, Q, &R, &Pcopy, A);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);

    } while (i--);
}

void exp_by_squaring_(fp* x, fp* y, uint64_t exp)
{
	fp result1, result2;
	fp_set(&result1, 1);
	fp_set(&result2, 1);

    while (exp)
    {
        if (exp & 1){
          fp_mul2(&result1, x);
          fp_mul2(&result2, y);
	}
	
        fp_sq1(x);
	      fp_sq1(y);
        exp >>= 1;
    }

    fp_cswap(&result1, x, 1);
    fp_cswap(&result2, y, 1);

}



void biquad_both_opt(fp *out,fp *outinv, fp *coeff, fp *coeffQ, const proj *A)
{
  fp t0, t1, t2;
  // coeffQ[0]=c(Qx+Qz)
  // coeffQ[1]=c(Qx-Qz)
  // coeffQ[2]= 4cQxQz
 // coeffQ[3]=2c(Qx^2+Qz^2)

  fp_mul3(&t0, &coeff[3], &coeffQ[1]); // t0 = (Px+Pz)*(Qx-Qz)
  fp_mul3(&t1, &coeff[4], &coeffQ[0]); // t1 = (Px-Pz)*(Qx+Qz)

  fp_add3(&out[0], &t0, &t1); // t2 = 2(PxQx-PzQz)
  fp_sq1(&out[0]); //4(PxQx-PzQz)^2
  //fp_mul2(&out[0],&A->z); //4C(PxQx-PzQz)^2
  

  fp_sub3(&out[2], &t1, &t0); // 2(PxQz-PzQx)
  fp_sq1(&out[2]); // 4(PxQz-PzQx)^2
  //fp_mul2(&out[2],&A->z);

  fp_mul3(&t0, &coeff[1], &coeffQ[3]); //t0= 4CPxPz*2c(Qx^2+Qz^2)
  fp_mul3(&t1, &coeff[0], &coeffQ[2]); // t1 = 4cQxQz*2C(Px^2+Pz^2)
  fp_mul3(&t2, &coeff[2], &coeffQ[2]); // t2 = 4APxPz* 4cQxQz

  fp_add3(&t0,&t0,&t1);
  fp_add3(&t0,&t0,&t2);
  fp_sub3(&out[1],&fp_0, &t0);

 
  outinv[1] = out[1];
  outinv[2] = out[0];
  outinv[0] = out[2];
}

void biquad_both_min(fp *out,fp *outinv,const proj *P,const proj *Q /*point*/,const proj *A)
{
  fp PxQx; fp_mul3(&PxQx,&P->x,&Q->x); //PxQx=PxQx
  fp PxQz; fp_mul3(&PxQz,&P->x,&Q->z); 
  fp PzQx; fp_mul3(&PzQx,&P->z,&Q->x);
  fp PzQz; fp_mul3(&PzQz,&P->z,&Q->z);
  fp PPQQ; fp_mul3(&PPQQ,&PxQx,&PzQz); //PPQQ=PxQxPzQz
  fp_add2(&PPQQ,&PPQQ); // PPQQ=2PxQxPzQz
  fp_mul2(&PPQQ,&A->x); // PPQQ=2APxQxPzQz

  fp s,t;

  fp_sub3(&s,&PxQx,&PzQz); // S = PxQx+PzQz
  fp_add3(&t,&PzQx,&PxQz); //t = PzQx+PxQz
  fp_mul3(&out[1],&s,&t); // out1 (PxQx+PzQz)(PzQx+PxQz)
  fp_mul2(&out[1],&A->z); // out1 C(PxQx+PzQz)(PzQx+PxQz)
  fp_add2(&out[1],&PPQQ);
  fp_add2(&out[1],&out[1]);
  fp_sub3(&out[1], &fp_0, &out[1]); /* XXX: push through other computations? */

  fp_sub3(&out[2],&PxQz,&PzQx);
  fp_sq1(&out[2]);
  fp_mul2(&out[2],&A->z);

  fp_add3(&out[0],&PxQx,&PzQz);
  fp_sq1(&out[0]);
  fp_mul2(&out[0],&A->z);

  /* ----- */
/*
  fp_add3(&s,&PxQz,&PzQx);
  fp_add3(&t,&PzQz,&PxQx);
  fp_mul3(&outinv[1],&s,&t);
  fp_mul2(&outinv[1],&A->z);
  fp_add2(&outinv[1],&PPQQ);
  fp_add2(&outinv[1],&outinv[1]);
  fp_neg1(&outinv[1]); XXX: push through other computations? */
  outinv[1] = out[1];
  fp_sub3(&outinv[1], &fp_0, &outinv[1]);
  outinv[2] = out[0];
  outinv[0] = out[2];
}


void biquad_both_min_opt(fp *out,fp *outinv, fp *out2,fp *outinv2, const proj *P, fp *coeff)
{
  // coeff[0]=2ACQxQz
  // coeff[1]=2ACQ2xQ2z
  // coeff[2] =CQx
  // coeff[3] =CQz
  // coeff[4] =CQ2x
  // coeff[5] =CQ2z
  fp PPQQ, PPQQ2;
  fp PxQx; fp_mul3(&PxQx,&P->x,&coeff[2]); //PxQx=CPxQx
  fp PxQz; fp_mul3(&PxQz,&P->x,&coeff[3]); 
  fp PzQx; fp_mul3(&PzQx,&P->z,&coeff[2]);
  fp PzQz; fp_mul3(&PzQz,&P->z,&coeff[3]);

  fp PxPz; fp_mul3(&PxPz,&P->x,&P->z); //PPQQ=PxQxPzQz
  fp_mul3(&PPQQ, &PxPz,&coeff[0]);  // 2ACPxQxPzQz

  fp PxQx2; fp_mul3(&PxQx2,&P->x,&coeff[4]); //PxQx=PxQx
  fp PxQz2; fp_mul3(&PxQz2,&P->x,&coeff[5]); 
  fp PzQx2; fp_mul3(&PzQx2,&P->z,&coeff[4]);
  fp PzQz2; fp_mul3(&PzQz2,&P->z,&coeff[5]);
  fp_mul3(&PPQQ2, &PxPz,&coeff[1]);  // 2ACPxQ2xPzQ2z


  fp s,t;
  fp_sub3(&out[2],&PxQz,&PzQx); // CPxQz-CPzQx
  fp_sq1(&out[2]); //s= c^2(PxQz-PzQx)^2 //out 2 = quad


  fp_add3(&out[0],&PxQx,&PzQz);
  fp_sq1(&out[0]); // constant

  fp_add3(&s, &PxQz,&PzQx);
  fp_sub3(&t, &PxQx,&PzQz);
  fp_mul3(&out[1],&s,&t); // out1 C^2(PxQx+PzQz)(PzQx+PxQz)
  fp_add2(&out[1],&PPQQ);
  fp_add2(&out[1],&out[1]);
  fp_sub3(&out[1], &fp_0, &out[1]); /* XXX: push through other computations? */

  outinv[1] = out[1];
  fp_sub3(&outinv[1], &fp_0, &outinv[1]);
  outinv[2] = out[0];
  outinv[0] = out[2];

  fp_sub3(&out2[2],&PxQz2,&PzQx2); // CPxQz-CPzQx
  fp_sq1(&out2[2]); //s= c^2(PxQz-PzQx)^2

  fp_add3(&out2[0],&PxQx2,&PzQz2);
  fp_sq1(&out2[0]);

  fp_add3(&s, &PxQz2,&PzQx2);
  fp_sub3(&t, &PxQx2,&PzQz2);
  fp_mul3(&out2[1],&s,&t); // out1 C^2(PxQx+PzQz)(PzQx+PxQz)
  fp_add2(&out2[1],&PPQQ2);
  fp_add2(&out2[1],&out2[1]);
  fp_sub3(&out2[1], &fp_0, &out2[1]); /* XXX: push through other computations? */

  outinv2[1] = out2[1];
  fp_sub3(&outinv2[1], &fp_0, &outinv2[1]);
  outinv2[2] = out2[0];
  outinv2[0] = out2[2];

}


void biquad_pm1_opt(fp *coeff, fp *outplus,fp *outminus,const proj *P,const proj *A)
{
  fp Pplus; fp_add3(&coeff[3],&P->x,&P->z); //plus=x+z
  fp Pminus; fp_sub3(&coeff[4],&P->x,&P->z); //pminus=x-z
  fp_sq2(&Pplus, &coeff[3]); //pplus=(x+z)^2
  fp_sq2(&Pminus, &coeff[4]); //pminus=(x-z)^2

  fp_mul3(&outplus[0],&Pminus,&A->z); // outp =C(x-z)^2
  outplus[2] = outplus[0];
  fp_mul3(&outminus[0],&Pplus,&A->z); // outm=C(x+z)^2
  outminus[2] = outminus[0];

  fp_add3(&coeff[0],&outplus[0],&outminus[0] ); // coeff[0]=2C(X^2+Z^2
  fp_sub3(&coeff[1], &outminus[0],&outplus[0] ); //4CXZ

  fp u;
  fp_sub3(&coeff[2], &Pplus, &Pminus); // coeff[2]=4XZ
  fp_mul2(&coeff[2],&A->x); // coeff[2=4xzA
  fp_sub3(&u, &fp_0, &coeff[2]);

 // fp_sub3(&u,&Pminus,&Pplus); // u=-4xz
//  fp_mul2(&u,&A->x); // u=-4xzA

  fp t;
  fp_add3(&t,&outminus[0],&outminus[0]); // t=2C(x^2+2xz+z^2)
  fp_sub3(&outplus[1],&u,&t); //outplus1=-4Axz-2C(x^2+2xz+z^2)

  fp_add3(&t,&outplus[0],&outplus[0]);
  fp_sub3(&outminus[1],&t,&u);
}

/* computes the isogeny with kernel point K of order k */
/* returns the new curve coefficient A and the image of P */
/* (obviously) not constant time in k */
void xISOG(proj *A, proj *P, proj const *K, uint64_t k)
{
    assert (k >= 3);
    assert (k % 2 == 1);

    fp tmp0, tmp1;
    fp T[4] = {K->z, K->x, K->x, K->z};
    proj Q;

    fp_mul3(&Q.x,  &P->x, &K->x);
    fp_mul3(&tmp0, &P->z, &K->z);
    fp_sub2(&Q.x,  &tmp0);

    fp_mul3(&Q.z,  &P->x, &K->z);
    fp_mul3(&tmp0, &P->z, &K->x);
    fp_sub2(&Q.z,  &tmp0);

    proj M[3] = {*K};
    xDBL(&M[1], A, K);

    for (uint64_t i = 1; i < k / 2; ++i) {

        if (i >= 2)
            xADD(&M[i % 3], &M[(i - 1) % 3], K, &M[(i - 2) % 3]);

        fp_mul3(&tmp0, &M[i % 3].x, &T[0]);
        fp_mul3(&tmp1, &M[i % 3].z, &T[1]);
        fp_add3(&T[0], &tmp0, &tmp1);

        fp_mul2(&T[1], &M[i % 3].x);

        fp_mul3(&tmp0, &M[i % 3].z, &T[2]);
        fp_mul3(&tmp1, &M[i % 3].x, &T[3]);
        fp_add3(&T[2], &tmp0, &tmp1);

        fp_mul2(&T[3], &M[i % 3].z);


        fp_mul3(&tmp0, &P->x, &M[i % 3].x);
        fp_mul3(&tmp1, &P->z, &M[i % 3].z);
        fp_sub2(&tmp0, &tmp1);
        fp_mul2(&Q.x,  &tmp0);

        fp_mul3(&tmp0, &P->x, &M[i % 3].z);
        fp_mul3(&tmp1, &P->z, &M[i % 3].x);
        fp_sub2(&tmp0, &tmp1);
        fp_mul2(&Q.z,  &tmp0);
    }

    fp_mul2(&T[0], &T[1]);
    fp_add2(&T[0], &T[0]); /* multiplication by 2 */

    fp_sq1(&T[1]);

    fp_mul2(&T[2], &T[3]);
    fp_add2(&T[2], &T[2]); /* multiplication by 2 */

    fp_sq1(&T[3]);

    /* Ax := T[1] * T[3] * Ax - 3 * Az * (T[1] * T[2] - T[0] * T[3]) */
    fp_mul3(&tmp0, &T[1], &T[2]);
    fp_mul3(&tmp1, &T[0], &T[3]);
    fp_sub2(&tmp0, &tmp1);
    fp_mul2(&tmp0, &A->z);
    fp_add3(&tmp1, &tmp0, &tmp0); fp_add2(&tmp0, &tmp1); /* multiplication by 3 */

    fp_mul3(&tmp1, &T[1], &T[3]);
    fp_mul2(&tmp1, &A->x);

    fp_sub3(&A->x, &tmp1, &tmp0);

    /* Az := Az * T[3]^2 */
    fp_sq1(&T[3]);
    fp_mul2(&A->z, &T[3]);

    /* X := X * Xim^2, Z := Z * Zim^2 */
    fp_sq1(&Q.x);
    fp_sq1(&Q.z);
    fp_mul2(&P->x, &Q.x);
    fp_mul2(&P->z, &Q.z);
}

// Isogeny eval : Montgomery curve
// coefficient : Edwards
void xISOG_hy(proj *A, proj *P, proj const *K, uint64_t k)
{
    assert (k >= 3);
    assert (k % 2 == 1);


    fp tmp0, tmp1, tmp2, Psum, Pdif;
    proj Q, Aed, prod;

    fp_add3(&Aed.z, &A->z, &A->z);  //compute twisted Edwards curve coefficients // C'=2C
    fp_add3(&Aed.x, &A->x, &Aed.z); // D'=A+2C
    fp_sub3(&Aed.z, &A->x, &Aed.z); //C' = A+2C-2C A
   
    fp_add3(&Psum, &P->x, &P->z);   //precomputations
    fp_sub3(&Pdif, &P->x, &P->z);

    fp_sub3(&prod.x, &K->x, &K->z);
    fp_add3(&prod.z, &K->x, &K->z);
    
    fp_mul3(&tmp1, &prod.x, &Psum);
    fp_mul3(&tmp0, &prod.z, &Pdif);
    fp_add3(&Q.x, &tmp0, &tmp1);
    fp_sub3(&Q.z, &tmp0, &tmp1);

    proj M[3] = {*K};
    xDBL(&M[1], A, K);

    for (uint64_t i = 1; i < k / 2; ++i) {

        if (i >= 2)
           xADD(&M[i % 3], &M[(i - 1) % 3], K, &M[(i - 2) % 3]);

	fp_sub3(&tmp1, &M[i % 3].x, &M[i % 3].z);
    	fp_add3(&tmp0, &M[i % 3].x, &M[i % 3].z);
	fp_mul2(&prod.x, &tmp1);
        fp_mul2(&prod.z, &tmp0);
    	fp_mul2(&tmp1, &Psum);
    	fp_mul2(&tmp0, &Pdif);
    	fp_add3(&tmp2, &tmp0, &tmp1);
	fp_mul2(&Q.x, &tmp2);
    	fp_sub3(&tmp2, &tmp0, &tmp1);
	fp_mul2(&Q.z, &tmp2);

    }


    // point evaluation
    fp_sq1(&Q.x);
    fp_sq1(&Q.z);
    fp_mul2(&P->x, &Q.x);
    fp_mul2(&P->z, &Q.z);

    //compute Aed.x^k, Aed.z^k
    exp_by_squaring_(&Aed.x, &Aed.z, k);

    //compute prod.x^8, prod.z^8
    fp_sq1(&prod.x);
    fp_sq1(&prod.x);
    fp_sq1(&prod.x);
    fp_sq1(&prod.z);
    fp_sq1(&prod.z);
    fp_sq1(&prod.z);

    //compute image curve parameters
    fp_mul2(&Aed.z, &prod.x);
    fp_mul2(&Aed.x, &prod.z);

    //compute Montgomery params
    fp_add3(&A->x, &Aed.x, &Aed.z);
    fp_sub3(&A->z, &Aed.x, &Aed.z);
    fp_add2(&A->x, &A->x);

}



void xISOG_sqrt(proj *A, proj *P, proj const *K, long long k)
{
    assert (k >= 3);
    assert (k % 2 == 1);


    long long bs = 0;
    long long gs = 0;

    steps(&bs,&gs,k);
 

    fp tmp0, tmp1, tmp2;

    proj Aed;
    fp_add3(&Aed.z, &A->z, &A->z);  //compute twisted Edwards curve coefficients
    fp_add3(&Aed.x, &A->x, &Aed.z);
    fp_sub3(&Aed.z, &A->x, &Aed.z);
   
    fp Psum, Pdif;
    fp coeffQ[4];
    
    fp_add3(&Psum, &P->x, &P->z);   //precomputations //Psum=x+z
    fp_sub3(&Pdif, &P->x, &P->z);   // Pdif= x-z
    fp_add3(&coeffQ[0], &P->x, &P->z);   //precomputations //Psum=Qx+Qz
    fp_sub3(&coeffQ[1], &P->x, &P->z);   // Pdif= Qx-Qz
    fp_sq2(&tmp0, &coeffQ[0]); //(Qx+Qz)^2
    fp_sq2(&tmp1, &coeffQ[1]); //(Qx-Qz)^2
    fp_add3(&coeffQ[3], &tmp0, &tmp1); // 2(Qx^2+Qz^2)
    fp_sub3(&coeffQ[2], &tmp0, &tmp1); // 4QzQz

    //included
    fp_mul3(&coeffQ[0], &coeffQ[0], &A->z); // c(Qx+Qz)
    fp_mul3(&coeffQ[1], &coeffQ[1], &A->z); // c(Qx-Qz)
    fp_mul3(&coeffQ[2], &coeffQ[2], &A->z); // c(Qx+Qz)
    fp_mul3(&coeffQ[3], &coeffQ[3], &A->z); // c(Qx-Qz)


    int Minit[k];
    proj M[k]; /* XXX: use less space */

    for (long long s = 0;s < k;++s) Minit[s] = 0;

    M[1] = *K; Minit[1] = 1;
    xDBL(&M[2], A, K); Minit[2] = 1;

    
      for (long long s = 3;s < k;++s) { //order k isogeny, 3,4,5
        if (s & 1) {
          long long i = s/2;
          assert(s == 2*i+1);
          if (i < bs) 
          { // if s = odd
            if (s == 3) 
            {
              assert(Minit[1]);
              assert(Minit[2]);
              xADD(&M[s],&M[2],&M[1],&M[1]); //M[3]=M[1]+M[2] // 3K //kernel
              Minit[s] = 1;
              continue;
            }
            assert(Minit[s-2]);
            assert(Minit[s-4]);
            assert(Minit[2]);
            xADD(&M[s],&M[s-2],&M[2],&M[s-4]);
            Minit[s] = 1;
            continue;
          }
        } 
        else 
        { // if s=even
          long long i = s/2-1;
          assert(s == 2*i+2);
          if (i < (k-1)/2-2*bs*gs) {
            if (s == 4) {
              assert(Minit[2]);
              xDBL(&M[s],A,&M[2]);
              Minit[s] = 1;
              continue;
            }
            assert(Minit[s-2]);
            assert(Minit[s-4]);
            assert(Minit[2]);
            xADD(&M[s],&M[s-2],&M[2],&M[s-4]);
            Minit[s] = 1;
            continue;
          }
        }
  
        if (bs > 0) {
          if (s == 2*bs) {
            assert(Minit[bs-1]);
            assert(Minit[bs+1]);
            assert(Minit[2]);
            xADD(&M[s],&M[bs+1],&M[bs-1],&M[2]);
            Minit[s] = 1;
            continue;
          } else if (s == 4*bs) {
            assert(Minit[2*bs]);
            xDBL(&M[s],A,&M[2*bs]);
            Minit[s] = 1;
            continue;
          } else if (s == 6*bs) {
            assert(Minit[2*bs]);
            assert(Minit[4*bs]);
            xADD(&M[s],&M[4*bs],&M[2*bs],&M[2*bs]);
            Minit[s] = 1;
            continue;
          } else if (s%(4*bs) == 2*bs) {
            long long j = s/(4*bs);
            assert(s == 2*bs*(2*j+1));
            if (j < gs) {
              assert(Minit[s-4*bs]);
              assert(Minit[s-8*bs]);
              assert(Minit[4*bs]);
              xADD(&M[s],&M[s-4*bs],&M[4*bs],&M[s-8*bs]);
              Minit[s] = 1;
              continue;
            }
          }
        }
      }
    

    proj Abatch;
    Abatch.x = fp_1;
    Abatch.z = fp_1;
    proj Qbatch;
    Qbatch.x = fp_1;
    Qbatch.z = fp_1;

      long long TIlen = 2*bs+poly_tree1size(bs);
      fp TI[TIlen];
  
      for (long long i = 0;i < bs;++i) {
        assert(Minit[2*i+1]);
        fp_sub3(&TI[2*i],&fp_0, &M[2*i+1].x);
        TI[2*i+1] = M[2*i+1].z;
      }
  
      poly_tree1(TI+2*bs,TI,bs);
  
      fp TP[3*gs];
      fp TPinv[3*gs];
      fp T1[3*gs];
      fp Tminus1[3*gs];
      fp coeff[5];

      for (long long j = 0;j < gs;++j) {
        assert(Minit[2*bs*(2*j+1)]);
        biquad_pm1_opt(coeff, T1+3*j,Tminus1+3*j,&M[2*bs*(2*j+1)],A);
        biquad_both_opt(TP+3*j,TPinv+3*j, coeff, coeffQ, A);

      }
  
      poly_multiprod2(TP,gs);
      poly_multiprod2(TPinv,gs);
      poly_multiprod2_selfreciprocal(T1,gs); // t1=1=c
      poly_multiprod2_selfreciprocal(Tminus1,gs); //tminus1=-1=1/c

      long long precompsize = poly_multieval_precomputesize(bs,2*gs+1);
      fp precomp[precompsize];
      poly_multieval_precompute(precomp,bs,2*gs+1,TI,TI+2*bs);

      fp v[bs];

      poly_multieval_postcompute(v,bs,TP,2*gs+1,TI,TI+2*bs,precomp);
      Qbatch.z = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Qbatch.z,&v[i]);

      poly_multieval_postcompute(v,bs,TPinv,2*gs+1,TI,TI+2*bs,precomp);
      Qbatch.x = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Qbatch.x,&v[i]);

      poly_multieval_postcompute(v,bs,T1,2*gs+1,TI,TI+2*bs,precomp);
      Abatch.x = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Abatch.x,&v[i]);

      poly_multieval_postcompute(v,bs,Tminus1,2*gs+1,TI,TI+2*bs,precomp);
      Abatch.z = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Abatch.z,&v[i]);

      for (long long i = 0;i < (k-1)/2-2*bs*gs;++i) {
        assert(Minit[2*i+2]);
        fp_sub3(&tmp1, &M[2*i+2].x, &M[2*i+2].z);
        fp_add3(&tmp0, &M[2*i+2].x, &M[2*i+2].z);
        fp_mul2(&Abatch.x,&tmp1);
        fp_mul2(&Abatch.z,&tmp0);
        fp_mul2(&tmp1, &Psum);
        fp_mul2(&tmp0, &Pdif);
        fp_add3(&tmp2, &tmp0, &tmp1);
        fp_mul2(&Qbatch.x, &tmp2);
        fp_sub3(&tmp2, &tmp0, &tmp1);
        fp_mul2(&Qbatch.z, &tmp2);
      }
    

    // point evaluation
    fp_sq1(&Qbatch.x);
    fp_sq1(&Qbatch.z);
    fp_mul2(&P->x, &Qbatch.x);
    fp_mul2(&P->z, &Qbatch.z);

    //compute Aed.x^k, Aed.z^k
    exp_by_squaring_(&Aed.x, &Aed.z, k);

    //compute prod.x^8, prod.z^8
    fp_sq1(&Abatch.x);
    fp_sq1(&Abatch.x);
    fp_sq1(&Abatch.x);
    fp_sq1(&Abatch.z);
    fp_sq1(&Abatch.z);
    fp_sq1(&Abatch.z);

    //compute image curve parameters
    fp_mul2(&Aed.z, &Abatch.x);
    fp_mul2(&Aed.x, &Abatch.z);

    //compute Montgomery params
    fp_add3(&A->x, &Aed.x, &Aed.z);
    fp_sub3(&A->z, &Aed.x, &Aed.z);
    fp_add2(&A->x, &A->x);
}



void xISOG_min_torsion(proj *A, proj *P, proj *P2, proj const *K, uint64_t k)
{
    assert (k >= 3);
    assert (k % 2 == 1);

    fp tmp0, tmp1;
    proj Q, Q2;

    fp_mul3(&Q.x,  &P->x, &K->x);   // Qx=XX_i
    fp_mul3(&tmp0, &P->z, &K->z);   // tmo0 = ZZ_i
    fp_add2(&Q.x,  &tmp0);          // Qx= XXi+ZZi

    fp_mul3(&Q.z,  &P->x, &K->z);   // Qz=XZi
    fp_mul3(&tmp0, &P->z, &K->x);   // tmo0=ZXi
    fp_sub2(&Q.z,  &tmp0);          // XZi-ZX_i

    fp_mul3(&Q2.x,  &P2->x, &K->x);   // Qx=XX_i
    fp_mul3(&tmp0, &P2->z, &K->z);   // tmo0 = ZZ_i
    fp_add2(&Q2.x,  &tmp0);          // Qx= XXi+ZZi

    fp_mul3(&Q2.z,  &P2->x, &K->z);   // Qz=XZi
    fp_mul3(&tmp0, &P2->z, &K->x);   // tmo0=ZXi
    fp_sub2(&Q2.z,  &tmp0);          // XZi-ZX_i

    proj M[3] = {*K};
    xDBL_min(&M[1], A, K);              // M[1]=2X_i

    for (uint64_t i = 1; i < k / 2; ++i) {

        if (i >= 2)
            xADD_min(&M[i % 3], &M[(i - 1) % 3], K, &M[(i - 2) % 3]);




        // computes isogeny (done)
        fp_mul3(&tmp0, &P->x, &M[i % 3].x); //Xx
        fp_mul3(&tmp1, &P->z, &M[i % 3].z); //ZZ
        fp_add2(&tmp0, &tmp1); //XX+ZZ
        fp_mul2(&Q.x,  &tmp0);

        fp_mul3(&tmp0, &P->x, &M[i % 3].z);
        fp_mul3(&tmp1, &P->z, &M[i % 3].x);
        fp_sub2(&tmp0, &tmp1);
        fp_mul2(&Q.z,  &tmp0);

        fp_mul3(&tmp0, &P2->x, &M[i % 3].x); //Xx
        fp_mul3(&tmp1, &P2->z, &M[i % 3].z); //ZZ
        fp_add2(&tmp0, &tmp1); //XX+ZZ
        fp_mul2(&Q2.x,  &tmp0);

        fp_mul3(&tmp0, &P2->x, &M[i % 3].z);
        fp_mul3(&tmp1, &P2->z, &M[i % 3].x);
        fp_sub2(&tmp0, &tmp1);
        fp_mul2(&Q2.z,  &tmp0);


    }



    /* X := X * Xim^2, Z := Z * Zim^2 */
    fp_sq1(&Q.x);
    fp_sq1(&Q.z);
    fp_mul2(&P->x, &Q.x);
    fp_mul2(&P->z, &Q.z);

    fp_sq1(&Q2.x);
    fp_sq1(&Q2.z);
    fp_mul2(&P2->x, &Q2.x);
    fp_mul2(&P2->z, &Q2.z);

    fp_sub3(&tmp0, &P2->z, &P2->x);
    fp_add3(&tmp1, &P2->z, &P2->x);
    fp_mul3(&A->x, &tmp0, &tmp1); //Z^2-X^2
    fp_mul3(&A->z, &P2->x, &P2->z); //t2=XZ
}



void xISOG_two_min_affine(proj *A, int64_t e, int64_t sign)
{
    fp t0, t1, t2, t3;

    if(e==0)
    {
        return;
    }
    if(sign)
    {
        fp_sub3(&A->x, &fp_0, &A->x);
    }
    fp_sq2(&t0, &A->x);                 // t0 = A^2
    fp_add3(&t0, &t0, &fp_4);           // t0 = A^2+4
    fp_sub3(&t1, &t0, &fp_1);           // t1 = A^2+3
    fp_sqrt_test(&t0, &t0);                       // t0 = sqrt(A^2+4)
    fp_mul3(&t0, &t0, &A->x);           // t0 = Asqrt(A^2+4)
    fp_sub3(&t0, &t0, &t1);             // t0 = Asqrt(A^2+4)-(A^2+3)
    fp_add3(&A->x, &t0, &t0);




    for(int i=2; i<e+1; i++)
    {
        fp_sq2(&t0, &A->x);                 // t0 = A^2
        fp_sub3(&t0, &t0, &fp_4);           // t0 = A^2-4
        fp_sqrt_test(&t0, &t0);                       // t0 = sqrt(A^2-4)
        fp_sub3(&t0, &A->x, &t0);           // t0 = A-sqrt(a^2-4)
        fp_mul3(&t0, &t0, &A->x);           // t0 = A(A-sqrt(A^2-4)) 
        fp_sub3(&t0, &fp_3, &t0);           // t0 = 3-A(A-sqrt(A^2-4))   
        fp_add3(&A->x, &t0, &t0);

    }

    fp_sq2(&t0, &A->x);             // t0 A^2
    fp_sub3(&t0, &t0, &fp_4);       // t0 =A^2-4
    fp_sqrt_test(&t0, &t0);                   // t0 = sqrt(A^2-4)
    fp_add3(&t1, &t0, &t0);         // t1 = 2sqrt(A^2-4)
    fp_add3(&t2, &t0, &A->x);       // t2 = A+sqrt(A^2-4)
    fp_mul3(&t0, &t2, &t0);         // t0 = sqrt(A^2-4)*(A+sqrt(A^2-4))
    fp_add3(&t3, &t0, &t0);         // t3 = 2sqrt(A^2-4)*(A+sqrt(A^2-4))
    fp_sqrtinv_test(&t3);
    fp_add3(&A->x, &t2, &t1);       // A' = A+3sqrt(A^2-4) 
    fp_mul3(&A->x, &A->x, &t3);   


    if(sign==0)
    {
        fp_sub3(&A->x, &fp_0, &A->x);
    }




}


void xISOG_min_to_mont(proj *A)
{
    fp t0, t1;
    
    fp_sq2(&t0, &A->x);             // t0 = A^2
    fp_add3(&t0, &t0, &fp_4);        // t0 = A^2+4  
    fp_sqrtinv_test(&t0);// t0 = sqrt(A^2+4)
    fp_add3(&t1, &A->x, &A->x);     // t1 = 2A
    fp_sub3(&t1, &fp_0, &t1);         // t1 = -2A
    fp_mul3(&A->x, &t0, &t1);
    
}

void xISOG_mont_to_min(proj *A)
{
    fp t0, t1;
    
    fp_sq2(&t0, &A->x);             // t0 = A^2

    fp_sub3(&t0, &fp_4, &t0);        // t0 = 4-A^2

    fp_sqrtinv_test(&t0);                   // t0 = sqrt(4-A^2)

    fp_add3(&t1, &A->x, &A->x);     // t1 = 2A
    fp_sub3(&t1, &fp_0, &t1);         // t1 = -2A
    fp_mul3(&A->x, &t0, &t1);
}


void xISOG_three_radical(proj *A,  proj const *K, int64_t e, int64_t sign)
{
    fp t0, t1, t2, t3, t4, t5;
    fp M, N, C, B;

   // Montgomery to tate

   // 3iso
   fp_copy(&M, &K->x);


    for(int64_t i=0; i<e; i++)
    {
        
      fp_sq2(&t0, &M); // t0 =t^2
      fp_mul3(&t1, &t0, &M); // t1 = t^3
      fp_add3(&N, &t1, &M); // N = t^3+t
      fp_add3(&t2, &N, &N);   // t2 = 2t^3+2t
      fp_add3(&t2, &t2, &N);   // t2 = 3t^3+3t
      fp_triexp(&N); // N= alpha

      fp_sub3(&t2, &t2, &M);    // t2 = 3t^3+2t

      fp_add3(&t3, &M, &N);      // t3 = t+a
      fp_mul3(&t3, &t3, &N);    // t3= a(a+t)
      fp_mul3(&t3, &t3, &M);    // t3 = at(a+t)

      fp_add3(&t0, &t3, &t3); // t0 = 2at(a+t)
      fp_add3(&t0, &t0, &t3); // t0 =3at(a+t)
      fp_add3(&M, &t0, &N);     // M = 3at(a+t)+a
      fp_add3(&M, &M, &t2);




    }
    // isogeny with kernel (0,0)



    // to mont
    fp_sq2(&t0, &M); // t0= t^2
    fp_mul3(&t1, &t0, &M);  // t1 = t^3

    fp_add3(&t1, &t1, &t1);  // t1 = 2t^3
    fp_add3(&t1, &t1, &t1);   // t1 = 4t^3

    fp_sub3(&t2, &t0, &fp_2);   // t2 = t^2-2
    fp_mul3(&t0, &t2, &t0);  // t0 =t^2(t^2-2)
    fp_add3(&t2, &t0, &t0);   // t2 = 2t(t^3-2)
    fp_add3(&t2, &t2, &t0);   // t2 = 3t(t^3-2)
    fp_sub3(&t2, &fp_1, &t2);   // t2 = 1-3t(t^3-2)
    fp_inv_test(&t1); 
    fp_mul3(&A->x, &t1, &t2);


}


void xISOG_sqrt_min(proj *A, proj *P, proj *P2, proj const *K, long long k)
{
    assert (k >= 3);
    assert (k % 2 == 1);


    long long bs = 0;
    long long gs = 0;

    steps(&bs,&gs,k);


    fp tmp0, tmp1, tmp2;
    fp t0, t1, t2;

    int Minit[k];
    proj M[k]; /* XXX: use less space */

    for (long long s = 0;s < k;++s) Minit[s] = 0;

    M[1] = *K; Minit[1] = 1;
    xDBL_min(&M[2], A, K); Minit[2] = 1;


      for (long long s = 3;s < k;++s) { //order k isogeny, 3,4,5
        if (s & 1) {
          long long i = s/2;
          assert(s == 2*i+1);
          if (i < bs) 
          { // if s = odd
            if (s == 3) 
            {
              assert(Minit[1]);
              assert(Minit[2]);
              xADD_min(&M[s],&M[2],&M[1],&M[1]); //M[3]=M[1]+M[2] // 3K //kernel
              Minit[s] = 1;
              continue;
            }
            assert(Minit[s-2]);
            assert(Minit[s-4]);
            assert(Minit[2]);
            xADD_min(&M[s],&M[s-2],&M[2],&M[s-4]);
            Minit[s] = 1;
            continue;
          }
        } 
        else 
        { // if s=even
          long long i = s/2-1;
          assert(s == 2*i+2);
          if (i < (k-1)/2-2*bs*gs) {
            if (s == 4) {
              assert(Minit[2]);
              xDBL_min(&M[s],A,&M[2]);
              Minit[s] = 1;
              continue;
            }
            assert(Minit[s-2]);
            assert(Minit[s-4]);
            assert(Minit[2]);
            xADD_min(&M[s],&M[s-2],&M[2],&M[s-4]);
            Minit[s] = 1;
            continue;
          }
        }
  
        if (bs > 0) {
          if (s == 2*bs) {
            assert(Minit[bs-1]);
            assert(Minit[bs+1]);
            assert(Minit[2]);
            xADD_min(&M[s],&M[bs+1],&M[bs-1],&M[2]);
            Minit[s] = 1;
            continue;
          } else if (s == 4*bs) {
            assert(Minit[2*bs]);
            xDBL_min(&M[s],A,&M[2*bs]);
            Minit[s] = 1;
            continue;
          } else if (s == 6*bs) {
            assert(Minit[2*bs]);
            assert(Minit[4*bs]);
            xADD_min(&M[s],&M[4*bs],&M[2*bs],&M[2*bs]);
            Minit[s] = 1;
            continue;
          } else if (s%(4*bs) == 2*bs) {
            long long j = s/(4*bs);
            assert(s == 2*bs*(2*j+1));
            if (j < gs) {
              assert(Minit[s-4*bs]);
              assert(Minit[s-8*bs]);
              assert(Minit[4*bs]);
              xADD_min(&M[s],&M[s-4*bs],&M[4*bs],&M[s-8*bs]);
              Minit[s] = 1;
              continue;
            }
          }
        }
      }
    

    proj Abatch;
    Abatch.x = fp_1;
    Abatch.z = fp_1;
    proj Qbatch;
    Qbatch.x = fp_1;
    Qbatch.z = fp_1;

    
      long long TIlen = 2*bs+poly_tree1size(bs);
      fp TI[TIlen];
  
      for (long long i = 0;i < bs;++i) {
        assert(Minit[2*i+1]);
        //fp_sub3(&TI[2*i],&base_p, &M[2*i+1].x);
        fp_sub3(&TI[2*i],&fp_0, &M[2*i+1].x);
        TI[2*i+1] = M[2*i+1].z;
      }
  
      poly_tree1(TI+2*bs,TI,bs);
  
      fp TP[3*gs];
      fp TPinv[3*gs];
      fp TP2[3*gs];
      fp TP2inv[3*gs];

  
      for (long long j = 0;j < gs;++j) {
        assert(Minit[2*bs*(2*j+1)]);
        biquad_both_min(TP+3*j,TPinv+3*j,&M[2*bs*(2*j+1)],P,A);
        biquad_both_min(TP2+3*j,TP2inv+3*j,&M[2*bs*(2*j+1)],P2,A);

      }
  
      poly_multiprod2(TP,gs);
      poly_multiprod2(TPinv,gs);
      poly_multiprod2(TP2,gs);
      poly_multiprod2(TP2inv,gs);

      long long precompsize = poly_multieval_precomputesize(bs,2*gs+1);
      fp precomp[precompsize];
      poly_multieval_precompute(precomp,bs,2*gs+1,TI,TI+2*bs);

      fp v[bs];

      poly_multieval_postcompute(v,bs,TP,2*gs+1,TI,TI+2*bs,precomp);
      Qbatch.z = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Qbatch.z,&v[i]);

      poly_multieval_postcompute(v,bs,TPinv,2*gs+1,TI,TI+2*bs,precomp);
      Qbatch.x = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Qbatch.x,&v[i]);

      poly_multieval_postcompute(v,bs,TP2,2*gs+1,TI,TI+2*bs,precomp);
      Abatch.z = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Abatch.z,&v[i]);

      poly_multieval_postcompute(v,bs,TP2inv,2*gs+1,TI,TI+2*bs,precomp);
      Abatch.x = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Abatch.x,&v[i]);

      for (long long i = 0;i < (k-1)/2-2*bs*gs;++i) {
        assert(Minit[2*i+2]);

        fp_mul3(&tmp0, &P->x, &M[2*i+2].x); //XX
        fp_mul3(&tmp1, &P->z, &M[2*i+2].z); //ZZ
        fp_add2(&tmp0, &tmp1); //XX+ZZ
        fp_mul2(&Qbatch.x,  &tmp0);

        fp_mul3(&tmp0, &P->x, &M[2*i+2].z); //XZ
        fp_mul3(&tmp1, &P->z, &M[2*i+2].x); //ZX
        fp_sub2(&tmp0, &tmp1);
        fp_mul2(&Qbatch.z,  &tmp0);

        fp_mul3(&tmp0, &P2->x, &M[2*i+2].x);
        fp_mul3(&tmp1, &P2->z, &M[2*i+2].z);
        fp_add2(&tmp0, &tmp1);
        fp_mul2(&Abatch.x,  &tmp0);

        fp_mul3(&tmp0, &P2->x, &M[2*i+2].z);
        fp_mul3(&tmp1, &P2->z, &M[2*i+2].x);
        fp_sub2(&tmp0, &tmp1);
        fp_mul2(&Abatch.z,  &tmp0);


      }
    

    // point evaluation
    fp_sq1(&Qbatch.x);
    fp_sq1(&Qbatch.z);
    fp_mul2(&P->x, &Qbatch.x);
    fp_mul2(&P->z, &Qbatch.z);

    //point evaluation for 2torsion
    fp_sq1(&Abatch.x);
    fp_sq1(&Abatch.z);
    fp_mul2(&P2->x, &Abatch.x);
    fp_mul2(&P2->z, &Abatch.z);

    //fp_sq2(&t0, &P2->z); //Z^2
    //fp_sq2(&t1, &P2->x); //X^2
    //fp_mul3(&A->z, &P2->x, &P2->z); //t2=XZ
    //fp_sub3(&A->x, &t0, &t1);

    fp_sub3(&t0, &P2->z, &P2->x);
    fp_add3(&t1, &P2->z, &P2->x);
    fp_mul3(&A->x, &t0, &t1); //Z^2-X^2
    fp_mul3(&A->z, &P2->x, &P2->z); //t2=XZ


}


void xISOG_sqrt_min_opt(proj *A, proj *P, proj *P2, proj const *K, long long k)
{
    assert (k >= 3);
    assert (k % 2 == 1);

    long long bs = 0;
    long long gs = 0;

    steps(&bs,&gs,k);


    fp tmp0, tmp1, tmp2;
    fp t0, t1, t2;
    fp coeff[6];
  
    fp_mul3(&coeff[0], &P->x, &P->z); 
    fp_mul3(&coeff[0], &coeff[0], &A->x);
    fp_mul3(&coeff[0], &coeff[0], &A->z);  

    fp_mul3(&coeff[1], &P2->x, &P2->z); 
    fp_mul3(&coeff[1], &coeff[1], &A->x); 
    fp_mul3(&coeff[1], &coeff[1], &A->z); 

    fp_add3(&coeff[0], &coeff[0],&coeff[0]);
    fp_add3(&coeff[1], &coeff[1],&coeff[1]);

    fp_mul3(&coeff[2], &A->z, &P->x);
    fp_mul3(&coeff[3], &A->z, &P->z);
    fp_mul3(&coeff[4], &A->z, &P2->x);
    fp_mul3(&coeff[5], &A->z, &P2->z);


    int Minit[k];
    proj M[k]; /* XXX: use less space */

    for (long long s = 0;s < k;++s) Minit[s] = 0;

    M[1] = *K; Minit[1] = 1;
    xDBL_min(&M[2], A, K); Minit[2] = 1;

    
      for (long long s = 3;s < k;++s) { //order k isogeny, 3,4,5
        if (s & 1) {
          long long i = s/2;
          assert(s == 2*i+1);
          if (i < bs) 
          { // if s = odd
            if (s == 3) 
            {
              assert(Minit[1]);
              assert(Minit[2]);
              xADD_min(&M[s],&M[2],&M[1],&M[1]); //M[3]=M[1]+M[2] // 3K //kernel
              Minit[s] = 1;
              continue;
            }
            assert(Minit[s-2]);
            assert(Minit[s-4]);
            assert(Minit[2]);
            xADD_min(&M[s],&M[s-2],&M[2],&M[s-4]);
            Minit[s] = 1;
            continue;
          }
        } 
        else 
        { // if s=even
          long long i = s/2-1;
          assert(s == 2*i+2);
          if (i < (k-1)/2-2*bs*gs) {
            if (s == 4) {
              assert(Minit[2]);
              xDBL_min(&M[s],A,&M[2]);
              Minit[s] = 1;
              continue;
            }
            assert(Minit[s-2]);
            assert(Minit[s-4]);
            assert(Minit[2]);
            xADD_min(&M[s],&M[s-2],&M[2],&M[s-4]);
            Minit[s] = 1;
            continue;
          }
        }
  
        if (bs > 0) {
          if (s == 2*bs) {
            assert(Minit[bs-1]);
            assert(Minit[bs+1]);
            assert(Minit[2]);
            xADD_min(&M[s],&M[bs+1],&M[bs-1],&M[2]);
            Minit[s] = 1;
            continue;
          } else if (s == 4*bs) {
            assert(Minit[2*bs]);
            xDBL_min(&M[s],A,&M[2*bs]);
            Minit[s] = 1;
            continue;
          } else if (s == 6*bs) {
            assert(Minit[2*bs]);
            assert(Minit[4*bs]);
            xADD_min(&M[s],&M[4*bs],&M[2*bs],&M[2*bs]);
            Minit[s] = 1;
            continue;
          } else if (s%(4*bs) == 2*bs) {
            long long j = s/(4*bs);
            assert(s == 2*bs*(2*j+1));
            if (j < gs) {
              assert(Minit[s-4*bs]);
              assert(Minit[s-8*bs]);
              assert(Minit[4*bs]);
              xADD_min(&M[s],&M[s-4*bs],&M[4*bs],&M[s-8*bs]);
              Minit[s] = 1;
              continue;
            }
          }
        }
      }
    
    proj Abatch;
    Abatch.x = fp_1;
    Abatch.z = fp_1;
    proj Qbatch;
    Qbatch.x = fp_1;
    Qbatch.z = fp_1;

    
      long long TIlen = 2*bs+poly_tree1size(bs);
      fp TI[TIlen];
  
      for (long long i = 0;i < bs;++i) {
        assert(Minit[2*i+1]);
        //fp_sub3(&TI[2*i],&base_p, &M[2*i+1].x);
        fp_sub3(&TI[2*i],&fp_0, &M[2*i+1].x);
        TI[2*i+1] = M[2*i+1].z;
      }
  
      poly_tree1(TI+2*bs,TI,bs);
  
      fp TP[3*gs];
      fp TPinv[3*gs];
      fp TP2[3*gs];
      fp TP2inv[3*gs];
  
  
      for (long long j = 0;j < gs;++j) {
        assert(Minit[2*bs*(2*j+1)]);
        biquad_both_min_opt(TP+3*j,TPinv+3*j,TP2+3*j,TP2inv+3*j,&M[2*bs*(2*j+1)],coeff);

      }
  
      poly_multiprod2(TP,gs);
      poly_multiprod2(TPinv,gs);
      poly_multiprod2(TP2,gs);
      poly_multiprod2(TP2inv,gs);

      long long precompsize = poly_multieval_precomputesize(bs,2*gs+1);
      fp precomp[precompsize];
      poly_multieval_precompute(precomp,bs,2*gs+1,TI,TI+2*bs);

      fp v[bs];

      poly_multieval_postcompute(v,bs,TP,2*gs+1,TI,TI+2*bs,precomp);
      Qbatch.z = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Qbatch.z,&v[i]);

      poly_multieval_postcompute(v,bs,TPinv,2*gs+1,TI,TI+2*bs,precomp);
      Qbatch.x = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Qbatch.x,&v[i]);

      poly_multieval_postcompute(v,bs,TP2,2*gs+1,TI,TI+2*bs,precomp);
      Abatch.z = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Abatch.z,&v[i]);

      poly_multieval_postcompute(v,bs,TP2inv,2*gs+1,TI,TI+2*bs,precomp);
      Abatch.x = v[0];
      for (long long i = 1;i < bs;++i) fp_mul2(&Abatch.x,&v[i]);

      for (long long i = 0;i < (k-1)/2-2*bs*gs;++i) {
        assert(Minit[2*i+2]);

        fp_mul3(&tmp0, &P->x, &M[2*i+2].x); //XX
        fp_mul3(&tmp1, &P->z, &M[2*i+2].z); //ZZ
        fp_add2(&tmp0, &tmp1); //XX+ZZ
        fp_mul2(&Qbatch.x,  &tmp0);

        fp_mul3(&tmp0, &P->x, &M[2*i+2].z); //XZ
        fp_mul3(&tmp1, &P->z, &M[2*i+2].x); //ZX
        fp_sub2(&tmp0, &tmp1);
        fp_mul2(&Qbatch.z,  &tmp0);

        fp_mul3(&tmp0, &P2->x, &M[2*i+2].x);
        fp_mul3(&tmp1, &P2->z, &M[2*i+2].z);
        fp_add2(&tmp0, &tmp1);
        fp_mul2(&Abatch.x,  &tmp0);

        fp_mul3(&tmp0, &P2->x, &M[2*i+2].z);
        fp_mul3(&tmp1, &P2->z, &M[2*i+2].x);
        fp_sub2(&tmp0, &tmp1);
        fp_mul2(&Abatch.z,  &tmp0);


      }
    


    // point evaluation
    fp_sq1(&Qbatch.x);
    fp_sq1(&Qbatch.z);
    fp_mul2(&P->x, &Qbatch.x);
    fp_mul2(&P->z, &Qbatch.z);

    //point evaluation for 2torsion
    fp_sq1(&Abatch.x);
    fp_sq1(&Abatch.z);
    fp_mul2(&P2->x, &Abatch.x);
    fp_mul2(&P2->z, &Abatch.z);


    fp_sub3(&t0, &P2->z, &P2->x);
    fp_add3(&t1, &P2->z, &P2->x);
    fp_mul3(&A->x, &t0, &t1); //Z^2-X^2
    fp_mul3(&A->z, &P2->x, &P2->z); //t2=XZ


}

