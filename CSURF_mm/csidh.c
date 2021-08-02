
#include <string.h>
#include <assert.h>

#include "uint.h"
#include "fp.h"
#include "mont.h"
#include "swindows.h"
#include "csidh.h"
#include "rng.h"

const public_key base = {0}; /* A = 0 */

/// Max exponent for CSIDH/Onuki's method 
uint8_t max_exp_csidh[NUM_PRIMES]={
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, //32
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, //64
    5, 5, 5, 5, 5, 5, 5, 5, 5,
};
#define const_CSIDH 365

//252.58 _ 342
#define const_CSURF 342
uint8_t max_exp_csurf[NUM_PRIMES+1]={
    3, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 
    4, 4, 4, 4, 4, 4, 4, 3, 3, 58
};


uint8_t max_exp_csurf_ori[NUM_PRIMES+1]={
    4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, //29
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
    4, 4, 4, 4, 4, 4, 4, 4, 4, 137
};

#define const_CRAD 335
uint8_t max_exp_crad_3[NUM_PRIMES+1]={
 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
    4, 4, 4, 4, 4, 3, 2, 2, 
     20, 30
};


/* TODO allow different encodings depending on parameters */
/* TODO waste less randomness */
void csidh_private(private_key *priv)
{
    memset(&priv->e, 0, sizeof(priv->e));
    for (size_t i = 0; i < NUM_PRIMES; ) {
        int8_t buf[64];
        randombytes(buf, sizeof(buf));
        for (size_t j = 0; j < sizeof(buf); ++j) {
            if (buf[j] <= MAX_EXPONENT && buf[j] >= -MAX_EXPONENT) {
                priv->e[i / 2] |= (buf[j] & 0xf) << i % 2 * 4;
                if (++i >= NUM_PRIMES)
                    break;
            }
        }
    }
}


void csidh_private_mod(uint8_t (*e)[NUM_PRIMES+1])
{
    int8_t buf[NUM_PRIMES+1];
    uint8_t tmp, sgn;

    randombytes(buf, sizeof(buf));
   for (size_t i = 0; i < NUM_PRIMES; i++ ) 
    {
       sgn=  (buf[i]&0x80)>>7;
       tmp = (buf[i]&0x7f)%(max_exp_csidh[i]+1);

       if(tmp==0)
       {
            e[0][i] = 0;
            e[1][i] = 0;     
       }
       else
       {
            e[sgn][i] = tmp;
            e[sgn^1][i] = 0;    
       }
       
    }

}


void csurf_private_original(uint8_t (*e)[NUM_PRIMES+1])
{
    int8_t buf[(NUM_PRIMES+1)<<1];
    uint8_t tmp, sgn;

    randombytes(buf, sizeof(buf));
   for (size_t i = 0; i < NUM_PRIMES+1; i++ ) 
    {
       sgn=  (buf[(i<<1)]&0x80)>>7;
       tmp = (buf[(i<<1)+1]&0xff)%(max_exp_csurf_ori[i]+1);

       if(tmp==0)
       {
            e[0][i] = 0;
            e[1][i] = 0;     
       }
       else
       {
            e[sgn][i] = tmp;
            e[sgn^1][i] = 0;    
       }
       
    }

}




void csurf_private_rad2(uint8_t (*e)[NUM_PRIMES+1])
{
    int8_t buf[(NUM_PRIMES+1)<<1];
    uint8_t tmp, sgn;

    randombytes(buf, sizeof(buf));
   for (size_t i = 0; i < NUM_PRIMES+1; i++ ) 
    {
       sgn=  (buf[(i<<1)]&0x80)>>7;
       tmp = (buf[(i<<1)+1]&0xff)%(max_exp_csurf[i]+1);

       if(tmp==0)
       {
            e[0][i] = 0;
            e[1][i] = 0;     
       }
       else
       {
            e[sgn][i] = tmp;
            e[sgn^1][i] = 0;    
       }
       
    }

}

void csurf_private_rad3(uint8_t (*e)[NUM_PRIMES+1])
{
    int8_t buf[(NUM_PRIMES+1)<<1];
    uint8_t tmp, sgn;

    randombytes(buf, sizeof(buf));
   for (size_t i = 0; i < NUM_PRIMES+1; i++ ) 
    {
       sgn=  (buf[(i<<1)]&0x80)>>7;
       tmp = (buf[(i<<1)+1]&0xff)%(max_exp_crad_3[i]+1);

       if(tmp==0)
       {
            e[0][i] = 0;
            e[1][i] = 0;     
       }
       else
       {
            e[sgn][i] = tmp;
            e[sgn^1][i] = 0;    
       }
       
    }

}



static bool validate_rec(proj *P, proj const *A, size_t lower, size_t upper, uint *order, bool *is_supersingular)
{
    assert(lower < upper);

    if (upper - lower == 1) {

        /* now P is [(p+1) / l_lower] times the original random point */
        /* we only gain information if this multiple is non-zero */

        if (memcmp(&P->z, &fp_0, sizeof(fp))) {

            uint tmp;
            uint_set(&tmp, primes[lower]);
            xMUL(P, A, P, &tmp);

            if (memcmp(&P->z, &fp_0, sizeof(fp))) {
                /* order does not divide p+1. */
                *is_supersingular = false;
                return true;
            }

            uint_mul3_64(order, order, primes[lower]);

            if (uint_sub3(&tmp, &four_sqrt_p, order)) { /* returns borrow */
                /* order > 4 sqrt(p), hence definitely supersingular */
                *is_supersingular = true;
                return true;
            }
        }

        /* inconclusive */
        return false;
    }

    size_t mid = lower + (upper - lower + 1) / 2;

    uint cl = uint_1, cu = uint_1;
    for (size_t i = lower; i < mid; ++i)
        uint_mul3_64(&cu, &cu, primes[i]);
    for (size_t i = mid; i < upper; ++i)
        uint_mul3_64(&cl, &cl, primes[i]);

    proj Q;

    xMUL(&Q, A, P, &cu);
    xMUL(P, A, P, &cl);

    /* start with the right half; bigger primes help more */
    return validate_rec(&Q, A, mid, upper, order, is_supersingular)
        || validate_rec(P, A, lower, mid, order, is_supersingular);
}

/* never accepts invalid keys. */
bool validate(public_key const *in)
{
    /* make sure the curve is nonsingular: A^2-4 != 0 */
    {
        uint dummy;
        if (!uint_sub3(&dummy, (uint *) &in->A, &p)) /* returns borrow */
            /* A >= p */
            return false;

        fp fp_pm2;
        fp_set(&fp_pm2, 2);
        if (!memcmp(&in->A, &fp_pm2, sizeof(fp)))
            /* A = 2 */
            return false;

        fp_sub3(&fp_pm2, &fp_0, &fp_pm2);
        if (!memcmp(&in->A, &fp_pm2, sizeof(fp)))
            /* A = -2 */
            return false;
    }

    const proj A = {in->A, fp_1};

    do {
        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        /* maximal 2-power in p+1 */
        xDBL(&P, &A, &P);
        xDBL(&P, &A, &P);

        bool is_supersingular;
        uint order = uint_1;

        if (validate_rec(&P, &A, 0, NUM_PRIMES, &order, &is_supersingular))
            return is_supersingular;

    /* P didn't have big enough order to prove supersingularity. */
    } while (1);
}

/* compute x^3 + Ax^2 + x */
static void montgomery_rhs(fp *rhs, fp const *A, fp const *x)
{
    fp tmp;
    *rhs = *x;
    fp_sq1(rhs);
    fp_mul3(&tmp, A, x);
    fp_add2(rhs, &tmp);
    fp_add2(rhs, &fp_1);
    fp_mul2(rhs, x);
}

static void montgomery_min_rhs(fp *rhs, fp const *A, fp const *x)
{
    fp tmp;
    *rhs = *x;
    fp_sq1(rhs); //x^2
    fp_mul3(&tmp, A, x); 
    fp_add2(rhs, &tmp); 
    fp_sub2(rhs, &fp_1); 
    fp_mul2(rhs, x); //
}


static void generate_two_torsion(proj *P1, proj *P2, fp const *A)
{
    fp rnd, rhs;
    bool done[2] = {false, false};
    bool sign; 
    do{
        fp_random(&rnd);
        montgomery_rhs(&rhs, A, &rnd);
        sign = !fp_issquare(&rhs); // 1=>squre 
        //if sign=0 -> square
        if(sign==0)
        {
            if(done[0]==true) continue;
            memcpy(&P1->x, &rnd, sizeof(fp));
            done[0] = true;
        }
        else
        {
            if(done[1]==true) continue;
            memcpy(&P2->x, &rnd, sizeof(fp));
            done[1] = true;        
        }
        


    }while (!(done[0] && done[1]));

}

static void generate_two_torsion_min(proj *P1, proj *P2, fp const *A)
{
    fp rnd, rhs;
    bool done[2] = {false, false};
    bool sign; 
    do{
        fp_random(&rnd);
        montgomery_min_rhs(&rhs, A, &rnd);
        sign = !fp_issquare(&rhs); // 1=>squre 
        //if sign=0 -> square
        if(sign==0)
        {
            if(done[0]==true) continue;
            memcpy(&P1->x, &rnd, sizeof(fp));
            done[0] = true;
        }
        else
        {
            if(done[1]==true) continue;
            memcpy(&P2->x, &rnd, sizeof(fp));
            done[1] = true;        
        }
        


    }while (!(done[0] && done[1]));

}



/* totally not constant-time. */
void action(public_key *out, public_key const *in,  uint8_t (*priv)[NUM_PRIMES+1])
{
    uint k[2];
    uint_set(&k[0], 1584); /* maximal 2-power in p+1 */
    uint_set(&k[1], 1584); /* maximal 2-power in p+1 */

    uint8_t e[2][NUM_PRIMES+1];
    memcpy(e[0], priv[0], sizeof(uint8_t)*(NUM_PRIMES+1));
    memcpy(e[1], priv[1], sizeof(uint8_t)*(NUM_PRIMES+1));

    for (size_t i = 0; i < NUM_PRIMES; i++) {

        if (e[1][i]!=0) {
            uint_mul3_64(&k[0], &k[0], primes[i]);

        }
        else if (e[0][i]!=0) {
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if(e[1][i]==0 && e[0][i]==0){
            uint_mul3_64(&k[0], &k[0], primes[i]);
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
    }

    proj A = {in->A, fp_1};

    bool done[2] = {false, false};

    do {

        assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        fp rhs;
        montgomery_rhs(&rhs, &A.x, &P.x);
        bool sign = !fp_issquare(&rhs);

        if (done[sign])
            continue;

        xMUL(&P, &A, &P, &k[sign]);

        done[sign] = true;

        for (size_t i = NUM_PRIMES - 1; i < NUM_PRIMES; --i) {

            if (e[sign][i]) {

                uint cof = uint_1;
                for (size_t j = 0; j < i; ++j)
                    if (e[sign][j])
                        uint_mul3_64(&cof, &cof, primes[j]);

                proj K;
                xMUL(&K, &A, &P, &cof);

                if (memcmp(&K.z, &fp_0, sizeof(fp))) {

                    xISOG(&A, &P, &K, primes[i]);

                    if (!--e[sign][i])
                        uint_mul3_64(&k[sign], &k[sign], primes[i]);

                }

            }

            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;

    } while (!(done[0] && done[1]));

    out->A = A.x;
}


void action_hy(public_key *out, public_key const *in, uint8_t (*priv)[NUM_PRIMES+1])
{
    uint k[2];
    uint_set(&k[0], 1584); /* maximal 2-power in p+1 */
    uint_set(&k[1], 1584); /* maximal 2-power in p+1 */

    uint8_t e[2][NUM_PRIMES+1];
   memcpy(e[0], priv[0], sizeof(uint8_t)*(NUM_PRIMES+1));
   memcpy(e[1], priv[1], sizeof(uint8_t)*(NUM_PRIMES+1));
    for (size_t i = 0; i < NUM_PRIMES; i++) {

        if (e[1][i]!=0) {
            uint_mul3_64(&k[0], &k[0], primes[i]);

        }
        else if (e[0][i]!=0) {
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if(e[1][i]==0 && e[0][i]==0){
            uint_mul3_64(&k[0], &k[0], primes[i]);
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
    }

    proj A = {in->A, fp_1};

    bool done[2] = {false, false};

    do {

        assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

        proj P;
        fp_random(&P.x);
        P.z = fp_1;
        fp rhs;
        montgomery_rhs(&rhs, &A.x, &P.x);
        bool sign = !fp_issquare(&rhs);

        if (done[sign])
            continue;
//ksign
        xMUL(&P, &A, &P, &k[sign]);
        done[sign] = true;

        for (size_t i = NUM_PRIMES-1; i < NUM_PRIMES; --i) {  //changed loop direction

            if (e[sign][i]) {

                uint cof = uint_1;
                for (size_t j = i - 1; j < NUM_PRIMES; --j)   //changed loop direction
                    if (e[sign][j])
                        uint_mul3_64(&cof, &cof, primes[j]);

                proj K;
               xMUL(&K, &A, &P, &cof);
              
                if (memcmp(&K.z, &fp_0, sizeof(fp))) {

                   xISOG_hy(&A, &P, &K, primes[i]);
   

                    if (!--e[sign][i])
                        uint_mul3_64(&k[sign], &k[sign], primes[i]);

                }

            }

            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;

    } while (!(done[0] && done[1]));

    out->A = A.x;
    
}

// Mont- curve using the 2torsion method and standard Velu
void action_min(public_key *out, public_key const *in, uint8_t (*priv)[NUM_PRIMES+1])
{
    uint k[2];
    uint_set(&k[0], 1584); /* maximal 2-power in p+1 */
    uint_set(&k[1], 1584); /* maximal 2-power in p+1 */

   uint8_t e[2][NUM_PRIMES+1];

   memcpy(e[0], priv[0], sizeof(uint8_t)*(NUM_PRIMES+1));
   memcpy(e[1], priv[1], sizeof(uint8_t)*(NUM_PRIMES+1));

    for (size_t i = 0; i < NUM_PRIMES; i++) {

        if (e[1][i]!=0) {
            uint_mul3_64(&k[0], &k[0], primes[i]);

        }
        else if (e[0][i]!=0) {
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if(e[1][i]==0 && e[0][i]==0){
            uint_mul3_64(&k[0], &k[0], primes[i]);
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
    }


    proj A = {in->A, fp_1};
    
    int64_t sign=0;
    if(e[1][NUM_PRIMES]!=0) sign = 1;

    xISOG_two_min_affine(&A, e[sign][NUM_PRIMES], sign);

        
    bool done[2] = {false, false};

    proj P2; //calc two torsion
    fp t0;
    fp_sq2(&t0, &A.x);             // t0 = A^2
    fp_add3(&t0, &t0, &fp_4);        // t0 = A^2+4
    fp_sqrt(&t0);                   // t0 = sqrt(A^2+4)
    fp_sub3(&P2.x, &t0, &A.x);
    P2.z=fp_2;

    do {

        assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        fp rhs;
        montgomery_min_rhs(&rhs, &A.x, &P.x);
        bool sign = !fp_issquare(&rhs);

        if (done[sign])
            continue;

        xMUL_min(&P, &A, &P, &k[sign]); // non_proj


        done[sign] = true;

        for (size_t i = NUM_PRIMES - 1; i < NUM_PRIMES; --i) {

            if (e[sign][i]) {

                uint cof = uint_1;
                for (size_t j = 0; j < i; ++j)
                    if (e[sign][j])
                        uint_mul3_64(&cof, &cof, primes[j]);

                proj K;
                xMUL_min(&K, &A, &P, &cof);

                if (memcmp(&K.z, &fp_0, sizeof(fp))) {

                    xISOG_min_torsion(&A, &P, &P2, &K, primes[i]);


                    if (!--e[sign][i])
                        uint_mul3_64(&k[sign], &k[sign], primes[i]);

                }

            }

            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;



    } while (!(done[0] && done[1]));

    out->A = A.x;

}


void action_hy_sqrt(public_key *out, public_key const *in, uint8_t (*priv)[NUM_PRIMES+1])
{
    uint k[2];
    uint_set(&k[0], 1584); /* maximal 2-power in p+1 */
    uint_set(&k[1], 1584); /* maximal 2-power in p+1 */



   uint8_t e[2][NUM_PRIMES+1];

   memcpy(e[0], priv[0], sizeof(uint8_t)*(NUM_PRIMES+1));
   memcpy(e[1], priv[1], sizeof(uint8_t)*(NUM_PRIMES+1));
    for (size_t i = 0; i < NUM_PRIMES; i++) {

        if (e[1][i]!=0) {
            uint_mul3_64(&k[0], &k[0], primes[i]);

        }
        else if (e[0][i]!=0) {
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if(e[1][i]==0 && e[0][i]==0){
            uint_mul3_64(&k[0], &k[0], primes[i]);
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
    }

    proj A = {in->A, fp_1};

    bool done[2] = {false, false};

    do {

        assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

        proj P;
        fp_random(&P.x);
        P.z = fp_1;
        fp rhs;
        montgomery_rhs(&rhs, &A.x, &P.x);
        bool sign = !fp_issquare(&rhs);

        if (done[sign])
            continue;
//ksign
        xMUL(&P, &A, &P, &k[sign]);
        done[sign] = true;

        for (size_t i = NUM_PRIMES-1; i < NUM_PRIMES; --i) {  //changed loop direction

            if (e[sign][i]) {

                uint cof = uint_1;
                for (size_t j = i - 1; j < NUM_PRIMES; --j)   //changed loop direction
                    if (e[sign][j])
                        uint_mul3_64(&cof, &cof, primes[j]);

                proj K;
               xMUL(&K, &A, &P, &cof);
              
                if (memcmp(&K.z, &fp_0, sizeof(fp))) {

                    if(primes[i]>100)
                    {
                        xISOG_sqrt(&A, &P, &K, primes[i]);
                    }                        
                    else
                    {
                        xISOG_hy(&A, &P, &K, primes[i]);
                    }
   

                    if (!--e[sign][i])
                        uint_mul3_64(&k[sign], &k[sign], primes[i]);

                }

            }

            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;

    } while (!(done[0] && done[1]));

    out->A = A.x;
    
}


void action_min_sqrt(public_key *out, public_key const *in, uint8_t (*priv)[NUM_PRIMES+1])
{
    uint k[2];
    uint_set(&k[0], 1584); /* maximal 2-power in p+1 */
    uint_set(&k[1], 1584); /* maximal 2-power in p+1 */

   uint8_t e[2][NUM_PRIMES+1];

   memcpy(e[0], priv[0], sizeof(uint8_t)*(NUM_PRIMES+1));
   memcpy(e[1], priv[1], sizeof(uint8_t)*(NUM_PRIMES+1));
    for (size_t i = 0; i < NUM_PRIMES; i++) {

        if (e[1][i]!=0) {
            uint_mul3_64(&k[0], &k[0], primes[i]);

        }
        else if (e[0][i]!=0) {
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if(e[1][i]==0 && e[0][i]==0){
            uint_mul3_64(&k[0], &k[0], primes[i]);
            uint_mul3_64(&k[1], &k[1], primes[i]);
        }
    }

    proj A = {in->A, fp_1};

    int64_t sign=0;
    if(e[1][NUM_PRIMES]!=0) sign = 1;


   xISOG_two_min_affine(&A, e[sign][NUM_PRIMES], sign);

        
    bool done[2] = {false, false};

    proj P2; //calc two torsion
    fp t0;
    fp_sq2(&t0, &A.x);             // t0 = A^2
    fp_add3(&t0, &t0, &fp_4);        // t0 = A^2+4
    fp_sqrt_test(&t0, &t0);                   // t0 = sqrt(A^2+4)
    fp_sub3(&P2.x, &t0, &A.x);
    P2.z=fp_2;

    do {

        //assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        fp rhs;
        montgomery_min_rhs(&rhs, &A.x, &P.x);
        bool sign = !fp_issquare(&rhs);

        if (done[sign])
            continue;

        xMUL_min(&P, &A, &P, &k[sign]); // non_proj

        done[sign] = true;

        for (size_t i = NUM_PRIMES - 1; i < NUM_PRIMES; --i) {

            if (e[sign][i]) {

                uint cof = uint_1;
                for (size_t j = 0; j < i; ++j)
                    if (e[sign][j])
                        uint_mul3_64(&cof, &cof, primes[j]);

                proj K;
                xMUL_min(&K, &A, &P, &cof); // A proj

                if (memcmp(&K.z, &fp_0, sizeof(fp))) {

                    if(primes[i]>60)
                    {
                        xISOG_sqrt_min_opt(&A, &P, &P2, &K, primes[i]);
                    }                        
                    else
                    {
                        xISOG_min_torsion(&A, &P, &P2, &K, primes[i]);
                    }
  

                    if (!--e[sign][i])
                        uint_mul3_64(&k[sign], &k[sign], primes[i]);

                }

            }

            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;

    } while (!(done[0] && done[1]));

    
    out->A = A.x;

}


void action_min_rad3(public_key *out, public_key const *in, uint8_t (*priv)[NUM_PRIMES+1])
{
    uint k[2];
    uint_set(&k[0], 4752); /* maximal 2-power in p+1 */
    uint_set(&k[1], 4752); /* maximal 2-power in p+1 */

   uint8_t e[2][NUM_PRIMES+1];

   memcpy(e[0], priv[0], sizeof(uint8_t)*(NUM_PRIMES+1));
   memcpy(e[1], priv[1], sizeof(uint8_t)*(NUM_PRIMES+1));
    for (size_t i = 0; i < NUM_PRIMES-1; i++) {

        if (e[1][i]!=0) {
            uint_mul3_64(&k[0], &k[0], primes_rad[i]);

        }
        else if (e[0][i]!=0) {
            uint_mul3_64(&k[1], &k[1], primes_rad[i]);
        }
        else if(e[1][i]==0 && e[0][i]==0){
            uint_mul3_64(&k[0], &k[0], primes_rad[i]);
            uint_mul3_64(&k[1], &k[1], primes_rad[i]);
        }
    }


    proj A = {in->A, fp_1};
    
    int64_t sign=0;
    if(e[1][NUM_PRIMES-1]!=0) sign = 1;

    xISOG_two_min_affine(&A, e[sign][NUM_PRIMES-1], sign);


    // select 3-torsion
   proj P, K;
    fp rhs;
    int found =0;
    bool sgn = false;

    if(e[1][NUM_PRIMES]!=0) sign = 1;
    else sign=0;

    // sgn = 0 if e <0
    // sgn = 1 if e > 0
    sgn = sign^1;
    while (found==0)
    {
        while (sgn == (sign^1))
        {
            fp_random(&P.x);
            P.z = fp_1;
            montgomery_min_rhs(&rhs, &A.x, &P.x);
            sgn = !fp_issquare(&rhs);
        }

        xMUL_min_affine(&K, &A, &P, &p_plus_1_t); // A proj
        if (memcmp(&K.z, &fp_0, sizeof(fp))) found =1;
        else
        {
            found =0;
            sgn = sign^1;
        }
        
    }




    fp_inv_test(&K.z);
    fp_mul2(&K.x, &K.z);
    K.z = fp_1;



    xISOG_three_radical(&A, &K, e[sign][NUM_PRIMES], sign);

        
    bool done[2] = {false, false};

    proj P2; //calc two torsion
    fp t0;
    fp_sq2(&t0, &A.x);             // t0 = A^2
    fp_add3(&t0, &t0, &fp_4);        // t0 = A^2+4
    fp_sqrt_test(&t0, &t0);                   // t0 = sqrt(A^2+4)
    fp_sub3(&P2.x, &t0, &A.x);
    P2.z=fp_2;

    do {

        assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        fp rhs;
        montgomery_min_rhs(&rhs, &A.x, &P.x);
        bool sign = !fp_issquare(&rhs);

        if (done[sign])
            continue;

        xMUL_min_affine(&P, &A, &P, &k[sign]); // non_proj

        done[sign] = true;

        for (size_t i = NUM_PRIMES - 2; i < NUM_PRIMES; --i){

            if (e[sign][i]) {

                uint cof = uint_1;
                for (size_t j = 0; j < i; ++j)
                    if (e[sign][j])
                        uint_mul3_64(&cof, &cof, primes_rad[j]);

                proj K;
                xMUL_min(&K, &A, &P, &cof); // A proj

                if (memcmp(&K.z, &fp_0, sizeof(fp))) {

  

                    if(primes_rad[i]>40)
                    {
                        xISOG_sqrt_min_opt(&A, &P, &P2, &K, primes_rad[i]);
                    }                        
                    else
                    {
                        xISOG_min_torsion(&A, &P, &P2, &K, primes_rad[i]);
                    }
  
       

                    if (!--e[sign][i])
                        uint_mul3_64(&k[sign], &k[sign], primes_rad[i]);

                }

            }

            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;

    } while (!(done[0] && done[1]));
 

    out->A = A.x;
}




// eval = Mont , coeff = Mont, original velu
bool csidh_original(public_key *out, public_key const *in, uint8_t (*priv)[NUM_PRIMES+1])
{
    /*
    if (!validate(in)) {
        fp_random(&out->A);
        return false;
    }
    */
    action(out, in, priv);
    return true;
}


// eval = Mont, coeff = Edwards, original velu
bool csidh_hybrid(public_key *out, public_key const *in, uint8_t (*priv)[NUM_PRIMES+1])
{
    /*
    if (!validate(in)) {
        fp_random(&out->A);
        return false;
    }
    */
    action_hy(out, in, priv);
    return true;
}


// eval = Mont, coeff = Edwards, sqrt_velu
bool csidh_hybrid_sqrt(public_key *out, public_key const *in, uint8_t (*priv)[NUM_PRIMES+1])
{
    /*
    if (!validate(in)) {
        fp_random(&out->A);
        return false;
    }
    */
    action_hy_sqrt(out, in, priv);
    return true;
}


// CSIDH on Montgomery minus curve, original velu
bool csidh_min(public_key *out, public_key const *in, uint8_t (*priv)[NUM_PRIMES+1])
{
    /*
    if (!validate(in)) {
        fp_random(&out->A);
        return false;
    }
    */
    action_min(out, in, priv);
    return true;
}

// CSIDH on Montgomery minus curve, sqrt velu
bool csidh_min_sqrt(public_key *out, public_key const *in, uint8_t (*priv)[NUM_PRIMES+1])
{
    /*
    if (!validate(in)) {
        fp_random(&out->A);
        return false;
    }
    */
    action_min_sqrt(out, in, priv);
    return true;
}

// CSIDH on Montgomery minus curve + radical 3 isog, sqrt velu
bool csidh_min_rad3_sqrt(public_key *out, public_key const *in, uint8_t (*priv)[NUM_PRIMES+1])
{
    /*
    if (!validate(in)) {
        fp_random(&out->A);
        return false;
    }
    */
    action_min_rad3(out, in, priv);
    return true;
}
