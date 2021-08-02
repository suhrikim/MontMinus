#ifndef MONT_H
#define MONT_H

#include "params.h"

void xDBL(proj *Q, proj const *A, proj const *P);
void xDBL_min(proj *Q, proj const *A, proj const *P);

void xADD(proj *S, proj const *P, proj const *Q, proj const *PQ);
void xADD_min(proj *S, proj const *P, proj const *Q, proj const *PQ);

void xDBLADD(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A);
void xDBLADD_min(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A);
void xDBLADD_min_affine(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A);


void xMUL(proj *Q, proj const *A, proj const *P, uint const *k);
void xMUL_min(proj *Q, proj const *A, proj const *P, uint const *k);
void xMUL_min_affine(proj *Q, proj const *A, proj const *P, uint const *k);


// original Velu : coeff : Mont, eval Mont
void xISOG(proj *A, proj *P, proj const *K, uint64_t k);

// original Velu : coeff : Edwards, eval Mont
void xISOG_hy(proj *A, proj *P, proj const *K, uint64_t k);

// sqrt Velu Mont
void xISOG_sqrt(proj *A, proj *P, proj const *K, long long k);

void xISOG_min_torsion(proj *A, proj *P, proj *P2, proj const *K, uint64_t k);

void xISOG_two_min_affine(proj *A, int64_t e, int64_t sign);
void xISOG_min_to_mont(proj *A);
void xISOG_mont_to_min(proj *A);
void xISOG_three_radical(proj *A,  proj const *K, int64_t e, int64_t sign);


void xISOG_sqrt_min_opt(proj *A, proj *P, proj *P2, proj const *K, long long k);
void xISOG_sqrt_min(proj *A, proj *P, proj *P2, proj const *K, long long k);

// Constant-time isogenies
void xISOG_opt_dummy(proj *A, proj *P1, proj *P2, proj const *K, uint64_t k, int mask);
void xISOG_const_mask(proj *A, proj *P1, proj *P2, proj const *K, uint64_t k, int mask);


void xISOG_two_min_const(proj *A, uint8_t e, uint8_t sign, uint8_t max);
void xISOG_sqrt_min_opt_const(proj *A, proj *TA, proj *P, proj *P2, proj const *K, uint64_t k, int mask);
void xISOG_const_mask_min(proj *A, proj *TA, proj *P1, proj *P2, proj const *K, uint64_t k, int mask);
void xISOG_three_radical_const(proj *A,  proj const *K, uint8_t e, uint8_t sign, uint8_t max);



#endif
