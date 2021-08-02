#ifndef PARAMS_H
#define PARAMS_H

#include <stdint.h>

#define LIMBS 8

typedef struct uint { uint64_t c[LIMBS]; } uint;
typedef struct fp { uint64_t c[LIMBS]; } fp;
typedef struct proj { struct fp x, z; } proj;

#define NUM_PRIMES 73
#define MAX_EXPONENT 5 /* (2*5+1)^74 is roughly 2^256 */


static const unsigned primes[NUM_PRIMES] = {
     3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
     67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137,  //31
     139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, // 47
     229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313,  //63
     317, 331, 337, 347, 349, 353, 359, 367, 373,
};


static const unsigned primes_rad[NUM_PRIMES-1] = {
        5,   7,   11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59, 
     61,  67,   71,  73,  79,  83,  89,  97, 101, 103, 107, 109, 113, 127, 131, 137,  //32
     139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, // 48
     229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313,  //64
     317, 331, 337, 347, 349, 353, 359, 367, 373,
};


#include "constants.h"

#endif
