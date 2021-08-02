#ifndef CSIDH_H
#define CSIDH_H

#include <stdbool.h>

#include "params.h"

typedef struct private_key {
    int8_t e[(NUM_PRIMES + 1) / 2]; /* packed int4_t */
} private_key;

typedef struct public_key {
    fp A; /* Montgomery coefficient: represents y^2 = x^3 + Ax^2 + x */
} public_key;

extern const public_key base;

void csidh_private(private_key *priv);
void csidh_private_mod(uint8_t (*e)[NUM_PRIMES+1]);
void csurf_private_rad2(uint8_t (*e)[NUM_PRIMES+1]);
void csurf_private_rad3(uint8_t (*e)[NUM_PRIMES+1]);
void csurf_private_original(uint8_t (*e)[NUM_PRIMES+1]);


bool csidh_original(public_key *out, public_key const *in, uint8_t (*priv)[NUM_PRIMES+1]);
bool csidh_hybrid(public_key *out, public_key const *in, uint8_t (*priv)[NUM_PRIMES+1]);
bool csidh_hybrid_sqrt(public_key *out, public_key const *in, uint8_t (*priv)[NUM_PRIMES+1]);
bool csidh_min(public_key *out, public_key const *in, uint8_t (*priv)[NUM_PRIMES+1]);
bool csidh_min_sqrt(public_key *out, public_key const *in, uint8_t (*priv)[NUM_PRIMES+1]);
bool csidh_min_rad3_sqrt(public_key *out, public_key const *in, uint8_t (*priv)[NUM_PRIMES+1]);



#endif
