
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <inttypes.h>

#include "rng.h"
#include "csidh.h"


#define OS_WIN       1
#define OS_LINUX     2

#if defined(__WINDOWS__)        // Microsoft Windows OS
    #define OS_TARGET OS_WIN
#else defined(__LINUX__)        // Linux OS
    #define OS_TARGET OS_LINUX 
#endif

int64_t cpucycles(void)
{ // Access system counter for benchmarking
#if (OS_TARGET == OS_WIN) && (TARGET == TARGET_AMD64 || TARGET == TARGET_x86)
    return __rdtsc();
#elif (OS_TARGET == OS_WIN) && (TARGET == TARGET_ARM)
    return __rdpmccntr64();
#elif (OS_TARGET == OS_LINUX) && (TARGET == TARGET_AMD64 || TARGET == TARGET_x86)
    unsigned int hi, lo;

    __asm__ __volatile__ ("rdtsc\n\t" : "=a" (lo), "=d"(hi));
    return ((int64_t)lo) | (((int64_t)hi) << 32);
#elif (OS_TARGET == OS_LINUX) && (TARGET == TARGET_ARM || TARGET == TARGET_ARM64)
    struct timespec time;

    clock_gettime(CLOCK_REALTIME, &time);
    return (int64_t)(time.tv_sec*1e9 + time.tv_nsec);
#else
    return 0;            
#endif
}

static __inline__ uint64_t rdtsc(void)
{
    uint32_t hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return lo | (uint64_t) hi << 32;
}

#define loop 10
/* defaults */

#ifndef BENCH_ITS
    #define BENCH_ITS 1000
#endif

#if !defined(BENCH_VAL) && !defined(BENCH_ACT)
    #define BENCH_VAL 1
    #define BENCH_ACT 1
#endif
#ifndef BENCH_VAL
    #define BENCH_VAL 0
#endif
#ifndef BENCH_ACT
    #define BENCH_ACT 0
#endif


const unsigned long its = BENCH_ITS;
const bool val = BENCH_VAL, act = BENCH_ACT;


const size_t stacksz = 0x8000;  /* 32k */

/* csidh.c */
bool validate(public_key const *in);
void action(public_key *out, public_key const *in, private_key const *priv);


void csidh_bench()
{

    unsigned long long  alice_keygen=0;
    int i;


    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;
    unsigned long long cycles1, cycles2;

    uint8_t Alice_e[2][NUM_PRIMES+1];
    uint8_t Bob_e[2][NUM_PRIMES+1];

    //public_key base = {fp_0}; /* A = 0 */

    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("BENCHMARKING CSIDH original \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    
  
    for(i=0; i<loop; i++)
    {
        csidh_private_mod(Alice_e);
        csidh_private_mod(Bob_e);

        cycles1 = cpucycles();
        csidh_original(&pub_alice, &base, Alice_e);
        csidh_original(&pub_bob, &base, Bob_e);
        csidh_original(&shared_alice, &pub_bob, Alice_e);
        csidh_original(&shared_bob, &pub_alice, Bob_e);
        cycles2 = cpucycles();
        alice_keygen = alice_keygen+((cycles2-cycles1)/4);
    }

    printf("  CSIDH group action runs in ............................... %10lld ", alice_keygen/loop);
    printf("\n");


    printf("\n\n");

}

void csidh_hy_bench()
{

    unsigned long long  alice_keygen=0;
    int i;


    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;
    unsigned long long cycles1, cycles2;

    uint8_t Alice_e[2][NUM_PRIMES+1];
    uint8_t Bob_e[2][NUM_PRIMES+1];

    //public_key base = {fp_0}; /* A = 0 */

    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("BENCHMARKING CSIDH hybrid original \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    
  
    for(i=0; i<loop; i++)
    {
        csidh_private_mod(Alice_e);
        csidh_private_mod(Bob_e);

        cycles1 = cpucycles();
        csidh_hybrid(&pub_alice, &base, Alice_e);
        csidh_hybrid(&pub_bob, &base, Bob_e);
        csidh_hybrid(&shared_alice, &pub_bob, Alice_e);
        csidh_hybrid(&shared_bob, &pub_alice, Bob_e);
        cycles2 = cpucycles();
        alice_keygen = alice_keygen+((cycles2-cycles1)/4);
    }

    printf("  CSIDH group action runs in ............................... %10lld ", alice_keygen/loop);
    printf("\n");


    printf("\n\n");

}

void csidh_min_original_bench()
{

    unsigned long long  alice_keygen=0;
    int i;


    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;
    unsigned long long cycles1, cycles2;

    uint8_t Alice_e[2][NUM_PRIMES+1];
    uint8_t Bob_e[2][NUM_PRIMES+1];

    //public_key base = {fp_0}; /* A = 0 */

    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("BENCHMARKING CSURF original parameter (original velu) \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    
  
    for(i=0; i<loop; i++)
    {
        csurf_private_original(Alice_e);
        csurf_private_original(Bob_e);

        cycles1 = cpucycles();
        csidh_min(&pub_alice, &base, Alice_e);
        csidh_min(&pub_bob, &base, Bob_e);
        csidh_min(&shared_alice, &pub_bob, Alice_e);
        csidh_min(&shared_bob, &pub_alice, Bob_e);
        cycles2 = cpucycles();
        alice_keygen = alice_keygen+((cycles2-cycles1)/4);
    }

    printf("  CSIDH group action runs in ............................... %10lld ", alice_keygen/loop);
    printf("\n");


    printf("\n\n");

}


void csidh_min_bench()
{

    unsigned long long  alice_keygen=0;
    int i;


    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;
    unsigned long long cycles1, cycles2;

    uint8_t Alice_e[2][NUM_PRIMES+1];
    uint8_t Bob_e[2][NUM_PRIMES+1];

    //public_key base = {fp_0}; /* A = 0 */

    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("BENCHMARKING CSURF original \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    
  
    for(i=0; i<loop; i++)
    {
        csurf_private_rad2(Alice_e);
        csurf_private_rad2(Bob_e);

        cycles1 = cpucycles();
        csidh_min(&pub_alice, &base, Alice_e);
        csidh_min(&pub_bob, &base, Bob_e);
        csidh_min(&shared_alice, &pub_bob, Alice_e);
        csidh_min(&shared_bob, &pub_alice, Bob_e);
        cycles2 = cpucycles();
        alice_keygen = alice_keygen+((cycles2-cycles1)/4);
    }

    printf("  CSIDH group action runs in ............................... %10lld ", alice_keygen/loop);
    printf("\n");


    printf("\n\n");

}


void csidh_hy_sqrt_bench()
{

    unsigned long long  alice_keygen=0;
    int i;


    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;
    unsigned long long cycles1, cycles2;

    uint8_t Alice_e[2][NUM_PRIMES+1];
    uint8_t Bob_e[2][NUM_PRIMES+1];

    //public_key base = {fp_0}; /* A = 0 */

    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("BENCHMARKING Hybrid-CSIDH using sqrt Velu \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    
  
    for(i=0; i<loop; i++)
    {
        csidh_private_mod(Alice_e);
        csidh_private_mod(Bob_e);

        cycles1 = cpucycles();
        csidh_hybrid_sqrt(&pub_alice, &base, Alice_e);
        csidh_hybrid_sqrt(&pub_bob, &base, Bob_e);
        csidh_hybrid_sqrt(&shared_alice, &pub_bob, Alice_e);
        csidh_hybrid_sqrt(&shared_bob, &pub_alice, Bob_e);
        cycles2 = cpucycles();
        alice_keygen = alice_keygen+((cycles2-cycles1)/4);
    }

    printf("  CSIDH group action runs in ............................... %10lld ", alice_keygen/loop);
    printf("\n");


    printf("\n\n");

}




void csidh_min_sqrt_bench()
{

    unsigned long long  alice_keygen=0;
    int i;


    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;
    unsigned long long cycles1, cycles2;

    uint8_t Alice_e[2][NUM_PRIMES+1];
    uint8_t Bob_e[2][NUM_PRIMES+1];

    //public_key base = {fp_0}; /* A = 0 */

    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("BENCHMARKING CSURF using sqrt Velu \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    
  
    for(i=0; i<loop; i++)
    {
        csurf_private_rad2(Alice_e);
        csurf_private_rad2(Bob_e);

        cycles1 = cpucycles();
        csidh_min_sqrt(&pub_alice, &base, Alice_e);
        csidh_min_sqrt(&pub_bob, &base, Bob_e);
        csidh_min_sqrt(&shared_alice, &pub_bob, Alice_e);
        csidh_min_sqrt(&shared_bob, &pub_alice, Bob_e);
        cycles2 = cpucycles();
        alice_keygen = alice_keygen+((cycles2-cycles1)/4);
    }

    printf("  CSIDH group action runs in ............................... %10lld ", alice_keygen/loop);
    printf("\n");


    printf("\n\n");

}

void csidh_rad3_sqrt_bench()
{

    unsigned long long  alice_keygen=0;
    int i;


    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;
    unsigned long long cycles1, cycles2;

    uint8_t Alice_e[2][NUM_PRIMES+1];
    uint8_t Bob_e[2][NUM_PRIMES+1];

    //public_key base = {fp_0}; /* A = 0 */

    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("BENCHMARKING CSURF using radical 3 isogeny and sqrt Velu \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    
  
    for(i=0; i<loop; i++)
    {
        csurf_private_rad3(Alice_e);
        csurf_private_rad3(Bob_e);

        cycles1 = cpucycles();
        csidh_min_rad3_sqrt(&pub_alice, &base, Alice_e);
        csidh_min_rad3_sqrt(&pub_bob, &base, Bob_e);
        csidh_min_rad3_sqrt(&shared_alice, &pub_bob, Alice_e);
        csidh_min_rad3_sqrt(&shared_bob, &pub_alice, Bob_e);
        cycles2 = cpucycles();
        alice_keygen = alice_keygen+((cycles2-cycles1)/4);
    }

    printf("  CSIDH group action runs in ............................... %10lld ", alice_keygen/loop);
    printf("\n");


    printf("\n\n");

}





void main ()
{
    // non-constant time


    csidh_bench();
    csidh_hy_bench();
    csidh_min_bench();
    csidh_min_original_bench();
 
    csidh_hy_sqrt_bench();
    csidh_min_sqrt_bench();
    csidh_rad3_sqrt_bench();
 




}