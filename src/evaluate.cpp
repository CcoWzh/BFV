#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <string.h>
#include "evaluate.h"
#include "unit.h"
#include "flint/fmpz_vec.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpz.h"

using namespace std;

void Encode(int m, fmpz_poly_t pm)
{
    int sgn;
    if (m < 0)
        sgn = -1;
    else
        sgn = 1;
    int num;
    m = abs(m);
    double bitnum = log2(m);
    int i;
    int bit;
    if (bitnum == (int)bitnum)
        num = (int)bitnum + 1;
    else
        num = (int)bitnum + 1;

    for (i = 0; i < num; i++)
    {
        bit = m % 2;
        fmpz_poly_set_coeff_ui(pm, i, bit);
        m = m / 2;
    }
    fmpz_poly_scalar_mul_si(pm, pm, sgn);
}

void SH_Encrypt(fmpz_poly_t m, fmpz_poly_t pk_p0, fmpz_poly_t pk_p1, fmpz_poly_t c0, fmpz_poly_t c1)
{
    fmpz_poly_t u;
    fmpz_poly_t e1;
    fmpz_poly_t e2;
    fmpz_poly_t temp;
    fmpz_poly_t temp1;
    fmpz_poly_init(u);
    fmpz_poly_init(e1);
    fmpz_poly_init(e2);
    fmpz_poly_init(temp);
    fmpz_poly_init(temp1);

    Gen_R2(u, 2);

    SampleD(e1, 2);
    SampleD(e2, 3);

    fmpz_poly_mul(temp, pk_p0, u);

    fmpz_poly_add(temp, temp, e1);

    fmpz_poly_scalar_mul_fmpz(temp1, m, Delta);

    fmpz_poly_add(temp, temp, temp1);

    fmpz_poly_scalar_smod_fmpz(temp, temp, q);

    fmpz_poly_set(c0, temp);

    fmpz_poly_mul(temp, pk_p1, u);

    fmpz_poly_add(temp, temp, e2);

    fmpz_poly_scalar_smod_fmpz(c1, temp, q);

    fmpz_poly_clear(u);
    fmpz_poly_clear(e1);
    fmpz_poly_clear(e2);
    fmpz_poly_clear(temp);
    fmpz_poly_clear(temp1);
}

void fmpz_poly_nearest_fmpz(fmpz_poly_t v)
{
    long deg = fmpz_poly_degree(v);
    int i;
    fmpz_t temp;
    fmpz_t temp1;
    fmpz_init(temp);
    fmpz_init(temp1);

    for (i = 0; i <= deg; i++)
    {
        fmpz_poly_get_coeff_fmpz(temp, v, i);

        fmpz_abs(temp1, temp);

        fmpz_fdiv_r(temp1, temp1, Delta);

        if (fmpz_equal(temp, Delta) == 1)
        {
            fmpz_set_ui(temp, 1);
        }
        else if (fmpz_equal(temp, neg_Delta) == 1)
        {
            fmpz_set_si(temp, -1);
        }
        else if (fmpz_is_zero(temp) == 1)
        {
            fmpz_set_ui(temp, 0);
        }

        else if (fmpz_sgn(temp) == -1 && fmpz_cmp(temp1, Delta_2) > 0)
        {
            fmpz_fdiv_q(temp, temp, Delta);
        }
        else if (fmpz_sgn(temp) == -1 && fmpz_cmp(temp1, Delta_2) <= 0)
        {
            fmpz_cdiv_q(temp, temp, Delta);
        }
        else if (fmpz_sgn(temp) == 1 && fmpz_cmp(temp1, Delta_2) >= 0)
        {
            fmpz_cdiv_q(temp, temp, Delta);
        }
        else
        {
            fmpz_fdiv_q(temp, temp, Delta);
        }

        fmpz_poly_set_coeff_fmpz(v, i, temp);
    }

    fmpz_clear(temp);
    fmpz_clear(temp1);
}

void SH_Decrypt(fmpz_poly_t sk, fmpz_poly_t c0, fmpz_poly_t c1,fmpz_poly_t dec_m)
{
    // fmpz_poly_t temp;
    fmpz_poly_init(dec_m);
    fmpz_poly_mul(dec_m, c1, sk);
    fmpz_poly_add(dec_m, dec_m, c0);
    fmpz_poly_scalar_smod_fmpz(dec_m, dec_m, q);
    fmpz_poly_nearest_fmpz(dec_m);
    fmpz_poly_scalar_smod_fmpz(dec_m, dec_m, t);

    // fmpz_poly_clear(dec_m);
}

// c0 = (c_00, c_01), c1 = (c_10, c_11), c = c0+c1
void SH_Add(fmpz_poly_t c_00, fmpz_poly_t c_01, fmpz_poly_t c_10, fmpz_poly_t c_11, fmpz_poly_t c0, fmpz_poly_t c1)
{
    fmpz_poly_t temp;
    fmpz_poly_init(temp);
    fmpz_poly_add(temp, c_00, c_10);
    fmpz_poly_scalar_smod_fmpz(c0, temp, q);

    fmpz_poly_add(temp, c_01, c_11);
    fmpz_poly_scalar_smod_fmpz(c1, temp, q);
    fmpz_poly_clear(temp);
}

void SH_Mul(fmpz_poly_t c_00, fmpz_poly_t c_01, fmpz_poly_t c_10, fmpz_poly_t c_11, fmpz_poly_t c0, fmpz_poly_t c1, fmpz_poly_t c2)
{
    fmpz_poly_t temp1;
    fmpz_poly_t temp2;

    fmpz_poly_init(temp1);
    fmpz_poly_init(temp2);

    fmpz_poly_mul(c0, c_00, c_10);
    fmpz_poly_nearest_fmpz(c0);

    fmpz_poly_mul(temp1, c_00, c_11);
    fmpz_poly_mul(temp2, c_01, c_10);
    fmpz_poly_add(c1, temp1, temp2);
    fmpz_poly_nearest_fmpz(c1);

    fmpz_poly_mul(c2, c_01, c_11);
    fmpz_poly_nearest_fmpz(c2);

    fmpz_poly_clear(temp1);
    fmpz_poly_clear(temp2);
}

void SH_DecMul(fmpz_poly_t sk, fmpz_poly_t c0, fmpz_poly_t c1, fmpz_poly_t c2)
{
    fmpz_poly_t temp;
    fmpz_poly_t temp2;
    fmpz_poly_init(temp);
    fmpz_poly_init(temp2);

    fmpz_poly_mul(temp, c1, sk);

    fmpz_poly_mul(temp2, c2, sk);
    fmpz_poly_mul(temp2, temp2, sk);

    fmpz_poly_add(temp, c0, temp);
    fmpz_poly_add(temp, temp, temp2);

    fmpz_poly_scalar_smod_fmpz(temp, temp, q);

    fmpz_poly_nearest_fmpz(temp);

    fmpz_poly_scalar_smod_fmpz(temp, temp, t);

    fmpz_poly_clear(temp);
    fmpz_poly_clear(temp2);
}

void Relinear(fmpz_poly_t *rlk_r0, fmpz_poly_t *rlk_r1, fmpz_poly_t c0, fmpz_poly_t c1, fmpz_poly_t c2, fmpz_poly_t _c0, fmpz_poly_t _c1)
{
    int i;
    int deg = fmpz_poly_degree(c2);

    fmpz_poly_t *c2_i;
    c2_i = (fmpz_poly_t *)malloc(sizeof(fmpz_poly_t) * l);

    fmpz_t temp;
    fmpz_t temp1;
    fmpz_poly_t temp2;
    fmpz_poly_t temp3;

    for (i = 0; i < l; i++)
        fmpz_poly_init(c2_i[i]);

    fmpz_init(temp);
    fmpz_init(temp1);
    fmpz_poly_init(temp2);
    fmpz_poly_init(temp3);

    for (i = 0; i <= deg; i++)
    {
        fmpz_poly_get_coeff_fmpz(temp, c2, i);

        fmpz_abs(temp1, temp);
        if (fmpz_cmp(temp1, T) < 0)
            fmpz_poly_set_coeff_fmpz(c2_i[0], i, temp);
        else
        {
            fmpz_fdiv_qr(temp, temp1, temp, T);
            fmpz_poly_set_coeff_fmpz(c2_i[0], i, temp1);
            fmpz_poly_set_coeff_fmpz(c2_i[1], i, temp);
        }
    }

    fmpz_poly_mul(temp2, rlk_r0[0], c2_i[0]);
    fmpz_poly_mul(temp3, rlk_r0[1], c2_i[1]);

    fmpz_poly_add(temp2, temp2, temp3);
    fmpz_poly_add(temp2, temp2, c0);
    fmpz_poly_scalar_smod_fmpz(_c0, temp2, q);

    fmpz_poly_mul(temp2, rlk_r1[0], c2_i[0]);
    fmpz_poly_mul(temp3, rlk_r1[1], c2_i[1]);

    fmpz_poly_add(temp2, temp2, temp3);
    fmpz_poly_add(temp2, temp2, c1);
    fmpz_poly_scalar_smod_fmpz(_c1, temp2, q);

    for (i = 0; i < l; i++)
        fmpz_poly_clear(c2_i[i]);

    fmpz_clear(temp);
    fmpz_clear(temp1);
    fmpz_poly_clear(temp2);
    fmpz_poly_clear(temp3);
}
