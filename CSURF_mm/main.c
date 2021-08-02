
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>

#include "fp.h"
#include "csidh.h"

void uint_print(uint const *x)
{
    for (size_t i = 8*LIMBS-1; i < 8*LIMBS; --i)
        printf("%02hhx", i[(unsigned char *) x->c]);
}

void fp_print(fp const *x)
{
    uint y;
    fp_dec(&y, x);
    uint_print(&y);
}

void test_csidh()
{
     clock_t t0, t1;

    public_key pub_alice, pub_bob;

    public_key shared_alice, shared_bob;
    uint8_t Alice_e[2][NUM_PRIMES+1];
    uint8_t Bob_e[2][NUM_PRIMES+1];

    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("TEST CSIDH Original \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    csidh_private_mod(Alice_e);
    csidh_private_mod(Bob_e); 


    t0 = clock();
    assert(csidh_original(&pub_alice, &base, Alice_e));
    t1 = clock();

    printf("Alice's public key    (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_original(&pub_bob, &base, Bob_e));
    t1 = clock();

    printf("Bob's public key      (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_bob.A);
    printf("\n\n");


    t0 = clock();
    assert(csidh_original(&shared_alice, &pub_bob, Alice_e));
    t1 = clock();

    printf("Alice's shared secret (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_original(&shared_bob, &pub_alice, Bob_e));
    t1 = clock();

    printf("Bob's shared secret   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_bob.A);
    printf("\n\n");

    printf("    ");
    if (memcmp(&shared_alice, &shared_bob, sizeof(public_key)))
        printf("\x1b[31mNOT EQUAL!\x1b[0m\n");
    else
        printf("\x1b[32mequal.\x1b[0m\n");
    printf("\n");

    printf("\n");
}


void test_csidh_hy()
{
     clock_t t0, t1;

    public_key pub_alice, pub_bob;

    public_key shared_alice, shared_bob;
    uint8_t Alice_e[2][NUM_PRIMES+1];
    uint8_t Bob_e[2][NUM_PRIMES+1];

    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("TEST CSIDH Hybrid \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    csidh_private_mod(Alice_e);
    csidh_private_mod(Bob_e); 


    t0 = clock();
    assert(csidh_hybrid(&pub_alice, &base, Alice_e));
    t1 = clock();

    printf("Alice's public key    (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_hybrid(&pub_bob, &base, Bob_e));
    t1 = clock();

    printf("Bob's public key      (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_bob.A);
    printf("\n\n");


    t0 = clock();
    assert(csidh_hybrid(&shared_alice, &pub_bob, Alice_e));
    t1 = clock();

    printf("Alice's shared secret (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_hybrid(&shared_bob, &pub_alice, Bob_e));
    t1 = clock();

    printf("Bob's shared secret   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_bob.A);
    printf("\n\n");

    printf("    ");
    if (memcmp(&shared_alice, &shared_bob, sizeof(public_key)))
        printf("\x1b[31mNOT EQUAL!\x1b[0m\n");
    else
        printf("\x1b[32mequal.\x1b[0m\n");
    printf("\n");

    printf("\n");
}



void test_csurf()
{
     clock_t t0, t1;

    public_key pub_alice, pub_bob;

    public_key shared_alice, shared_bob;
    uint8_t Alice_e[2][NUM_PRIMES+1];
    uint8_t Bob_e[2][NUM_PRIMES+1];

    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("TEST CSURF original \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    csurf_private_rad2(Alice_e);
    csurf_private_rad2(Bob_e); 


    t0 = clock();
    assert(csidh_min(&pub_alice, &base, Alice_e));
    t1 = clock();

    printf("Alice's public key    (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_min(&pub_bob, &base, Bob_e));
    t1 = clock();

    printf("Bob's public key      (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_bob.A);
    printf("\n\n");


    t0 = clock();
    assert(csidh_min(&shared_alice, &pub_bob, Alice_e));
    t1 = clock();

    printf("Alice's shared secret (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_min(&shared_bob, &pub_alice, Bob_e));
    t1 = clock();

    printf("Bob's shared secret   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_bob.A);
    printf("\n\n");

    printf("    ");
    if (memcmp(&shared_alice, &shared_bob, sizeof(public_key)))
        printf("\x1b[31mNOT EQUAL!\x1b[0m\n");
    else
        printf("\x1b[32mequal.\x1b[0m\n");
    printf("\n");

    printf("\n");
}


void test_hybrid_sqrt()
{
     clock_t t0, t1;

    public_key pub_alice, pub_bob;

    public_key shared_alice, shared_bob;
    uint8_t Alice_e[2][NUM_PRIMES+1];
    uint8_t Bob_e[2][NUM_PRIMES+1];

    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("TEST Hybrid-CSIDH using sqrt Velu \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    csidh_private_mod(Alice_e);
    csidh_private_mod(Bob_e); 


    t0 = clock();
    assert(csidh_hybrid_sqrt(&pub_alice, &base, Alice_e));
    t1 = clock();

    printf("Alice's public key    (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_hybrid_sqrt(&pub_bob, &base, Bob_e));
    t1 = clock();

    printf("Bob's public key      (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_bob.A);
    printf("\n\n");


    t0 = clock();
    assert(csidh_hybrid_sqrt(&shared_alice, &pub_bob, Alice_e));
    t1 = clock();

    printf("Alice's shared secret (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_hybrid_sqrt(&shared_bob, &pub_alice, Bob_e));
    t1 = clock();

    printf("Bob's shared secret   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_bob.A);
    printf("\n\n");

    printf("    ");
    if (memcmp(&shared_alice, &shared_bob, sizeof(public_key)))
        printf("\x1b[31mNOT EQUAL!\x1b[0m\n");
    else
        printf("\x1b[32mequal.\x1b[0m\n");
    printf("\n");

    printf("\n");
}


void test_csurf_sqrt()
{
     clock_t t0, t1;

    public_key pub_alice, pub_bob;

    public_key shared_alice, shared_bob;
    uint8_t Alice_e[2][NUM_PRIMES+1];
    uint8_t Bob_e[2][NUM_PRIMES+1];

    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("TEST CSURF using sqrt Velu \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    csurf_private_rad2(Alice_e);
    csurf_private_rad2(Bob_e); 


    t0 = clock();
    assert(csidh_min_sqrt(&pub_alice, &base, Alice_e));
    t1 = clock();

    printf("Alice's public key    (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_min_sqrt(&pub_bob, &base, Bob_e));
    t1 = clock();

    printf("Bob's public key      (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_bob.A);
    printf("\n\n");


    t0 = clock();
    assert(csidh_min_sqrt(&shared_alice, &pub_bob, Alice_e));
    t1 = clock();

    printf("Alice's shared secret (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_min_sqrt(&shared_bob, &pub_alice, Bob_e));
    t1 = clock();

    printf("Bob's shared secret   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_bob.A);
    printf("\n\n");

    printf("    ");
    if (memcmp(&shared_alice, &shared_bob, sizeof(public_key)))
        printf("\x1b[31mNOT EQUAL!\x1b[0m\n");
    else
        printf("\x1b[32mequal.\x1b[0m\n");
    printf("\n");

    printf("\n");
}

void test_csurf_rad3_sqrt()
{
     clock_t t0, t1;

    public_key pub_alice, pub_bob;

    public_key shared_alice, shared_bob;
    uint8_t Alice_e[2][NUM_PRIMES+1];
    uint8_t Bob_e[2][NUM_PRIMES+1];

    printf("\n");

    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 
    printf("TEST CSURF using radical 3 and sqrt Velu \n"); 
    printf("\n--------------------------------------------------------------------------------------------------------\n\n"); 


    csurf_private_rad3(Alice_e);
    csurf_private_rad3(Bob_e); 


    t0 = clock();
    assert(csidh_min_rad3_sqrt(&pub_alice, &base, Alice_e));
    t1 = clock();

    printf("Alice's public key    (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_min_rad3_sqrt(&pub_bob, &base, Bob_e));
    t1 = clock();

    printf("Bob's public key      (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&pub_bob.A);
    printf("\n\n");


    t0 = clock();
    assert(csidh_min_rad3_sqrt(&shared_alice, &pub_bob, Alice_e));
    t1 = clock();

    printf("Alice's shared secret (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_alice.A);
    printf("\n\n");

    t0 = clock();
    assert(csidh_min_rad3_sqrt(&shared_bob, &pub_alice, Bob_e));
    t1 = clock();

    printf("Bob's shared secret   (%7.3lf ms):\n  ", 1000. * (t1 - t0) / CLOCKS_PER_SEC);
    fp_print(&shared_bob.A);
    printf("\n\n");

    printf("    ");
    if (memcmp(&shared_alice, &shared_bob, sizeof(public_key)))
        printf("\x1b[31mNOT EQUAL!\x1b[0m\n");
    else
        printf("\x1b[32mequal.\x1b[0m\n");
    printf("\n");

    printf("\n");
}



void main()
{

    // test nonconst
    
    test_csidh();
    test_csidh_hy();
    test_csurf();
    test_hybrid_sqrt();
    test_csurf_sqrt();
    test_csurf_rad3_sqrt();
    

    //test const


}