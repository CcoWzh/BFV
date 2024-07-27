#include <stdlib.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "key.h"
#include "flint/fmpz_vec.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpz.h"

void secret_key_Gen(fmpz_poly_t sk)
{
   Gen_R2(sk, 1);
}

void public_key_Gen(fmpz_poly_t sk, fmpz_poly_t pk_p0, fmpz_poly_t pk_p1)
{
   // pk = ( [-(as+e)]_q, a ), a \in R_q, s \in R_2,  e \in D_{Z, sigma}^d
   fmpz_poly_t e;
   fmpz_poly_t temp;
   fmpz_poly_init(e);
   fmpz_poly_init(temp);

   SampleD(e, 1);
   Gen_Rq(pk_p1);

   fmpz_poly_mul(temp, pk_p1, sk);
   fmpz_poly_add(temp, temp, e);
   fmpz_poly_neg(temp, temp);
   fmpz_poly_scalar_smod_fmpz(temp, temp, q);
   fmpz_poly_set(pk_p0, temp);

//    fmpz_poly_print(pk_p1);

   fmpz_poly_clear(e);
   fmpz_poly_clear(temp);
}

void ReKey_Gen(fmpz_poly_t sk, fmpz_poly_t *rlk_r0, fmpz_poly_t *rlk_r1)
{
   fmpz_poly_t e;
   int i;
   fmpz_poly_t temp1;
   fmpz_t temp2;
   fmpz_poly_t temp3;
   fmpz_poly_init(e);
   fmpz_poly_init(temp1);
   fmpz_init(temp2);
   fmpz_poly_init(temp3);

   for (i = 0; i < l; i++)
   {
      Gen_Rq(rlk_r1[i]);
      SampleD(e, i + 9);
      fmpz_poly_mul(temp1, rlk_r1[i], sk);
      fmpz_poly_add(temp1, temp1, e);
      fmpz_poly_neg(temp1, temp1);
      fmpz_pow_ui(temp2, T, i);
      fmpz_poly_mul(temp3, sk, sk);
      fmpz_poly_scalar_mul_fmpz(temp3, temp3, temp2);
      fmpz_poly_add(temp1, temp1, temp3);
      fmpz_poly_scalar_smod_fmpz(rlk_r0[i], temp1, q);
   }

   fmpz_poly_clear(temp1);
   fmpz_clear(temp2);
   fmpz_poly_clear(temp3);
   fmpz_poly_clear(e);
}

// void SH_Keygen(fmpz_poly_t sk, fmpz_poly_t pk_p0, fmpz_poly_t pk_p1, fmpz_poly_t rlk_r0, fmpz_poly_t rlk_r1)
// {
//    secret_key_Gen(sk);               // 生成私钥
//    public_key_Gen(sk, pk_p0, pk_p1); // 生成公钥
//    ReKey_Gen(sk, rlk_r0, rlk_r1);    // 生成Relinear key
// }