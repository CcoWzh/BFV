#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <string.h>
#include "key.h"
#include "parameter.h"
#include "evaluate.h"
#include "unit.h"
#include "flint/fmpz_vec.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpz.h"

using namespace std;

int main()
{
   fmpz_poly_t pk_p0;
   fmpz_poly_t pk_p1;
   fmpz_poly_t sk;
   fmpz_poly_t *rlk_r0;
   fmpz_poly_t *rlk_r1;

   int i;
   fmpz_poly_t m0;
   fmpz_poly_t m1;
   fmpz_poly_t dec_m;
   fmpz_poly_t c_00;
   fmpz_poly_t c_01;
   fmpz_poly_t c_10;
   fmpz_poly_t c_11;
   fmpz_poly_t c0;
   fmpz_poly_t c1;
   fmpz_poly_t c2;
   fmpz_poly_t _c0;
   fmpz_poly_t _c1;
   struct timespec tpstart;
   struct timespec tpend;
   struct timespec tp;
   long timedif;

   fmpz_init(q);
   fmpz_init(neg_Delta);
   fmpz_init(Delta_2);
   fmpz_init(t);
   fmpz_init(Delta);
   fmpz_init(T);

   fmpz_poly_init(pk_p0);
   fmpz_poly_init(pk_p1);
   fmpz_poly_init(sk);
   fmpz_poly_init(m0);
   fmpz_poly_init(m1);
   fmpz_poly_init(dec_m);
   fmpz_poly_init(c_00);
   fmpz_poly_init(c_01);
   fmpz_poly_init(c_10);
   fmpz_poly_init(c_11);
   fmpz_poly_init(c0);
   fmpz_poly_init(c1);
   fmpz_poly_init(c2);
   fmpz_poly_init(_c0);
   fmpz_poly_init(_c1);

   fmpz_t result;
   fmpz_init(result);

   System_Param();

   rlk_r0 = (fmpz_poly_t *)malloc(sizeof(fmpz_poly_t) * (l));
   rlk_r1 = (fmpz_poly_t *)malloc(sizeof(fmpz_poly_t) * (l));

   for (i = 0; i < l; i++)
   {
      fmpz_poly_init(rlk_r0[i]);
      fmpz_poly_init(rlk_r1[i]);
   }

   cout << "generate sk pk rlk ..." << endl;
   secret_key_Gen(sk);
   public_key_Gen(sk, pk_p0, pk_p1);
   ReKey_Gen(sk, rlk_r0, rlk_r1);

   Encode(5, m0);
   Encode(8, m1);

   SH_Encrypt(m0, pk_p0, pk_p1, c_00, c_01);
   SH_Encrypt(m1, pk_p0, pk_p1, c_10, c_11);

   SH_Decrypt(sk, c_00, c_01, dec_m);
   poly_to_num(dec_m, result);
   cout << "decrypt m0 is ";
   fmpz_print(result);
   flint_printf("\n");

   SH_Decrypt(sk, c_10, c_11, dec_m);
   poly_to_num(dec_m, result);
   cout << "decrypt m1 is ";
   fmpz_print(result);
   flint_printf("\n");

   SH_Add(c_00, c_01, c_10, c_11, c0, c1);
   SH_Decrypt(sk, c0, c1, dec_m);
   poly_to_num(dec_m, result);
   cout << "decrypt m0 + m1 = ";
   fmpz_print(result);
   flint_printf("\n");

   SH_Mul(c_00, c_01, c_10, c_11, c0, c1, c2);
   // SH_DecMul(sk, c0, c1, c2);
   Relinear(rlk_r0, rlk_r1, c0, c1, c2, _c0, _c1);
   SH_Decrypt(sk, _c0, _c1, dec_m);
   poly_to_num(dec_m, result);
   cout << "decrypt m0 * m1 = ";
   fmpz_print(result);
   flint_printf("\n");

   fmpz_poly_clear(sk);
   fmpz_poly_clear(pk_p0);
   fmpz_poly_clear(pk_p1);
   fmpz_clear(q);
   fmpz_clear(neg_Delta);
   fmpz_clear(Delta_2);
   fmpz_clear(t);
   fmpz_clear(Delta);
   fmpz_clear(T);
   fmpz_clear(result);
   fmpz_poly_clear(c_00);
   fmpz_poly_clear(c_01);
   fmpz_poly_clear(c_10);
   fmpz_poly_clear(c_11);
   fmpz_poly_clear(c0);
   fmpz_poly_clear(c1);
   fmpz_poly_clear(c2);
   fmpz_poly_clear(_c0);
   fmpz_poly_clear(_c1);
   fmpz_poly_clear(m0);
   fmpz_poly_clear(m1);
   fmpz_poly_clear(dec_m);

   for (i = 0; i < l; i++)
   {
      fmpz_poly_clear(rlk_r0[i]);
      fmpz_poly_clear(rlk_r1[i]);
   }

   return 0;
}