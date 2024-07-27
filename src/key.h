#pragma once
#include "parameter.h"
#include "flint/fmpz_vec.h"
#include "flint/fmpz_poly.h"

void secret_key_Gen(fmpz_poly_t sk);
void public_key_Gen(fmpz_poly_t sk, fmpz_poly_t pk_p0, fmpz_poly_t pk_p1);
void ReKey_Gen(fmpz_poly_t sk, fmpz_poly_t *rlk_r0, fmpz_poly_t *rlk_r1);

// void SH_Keygen(fmpz_poly_t sk, fmpz_poly_t pk_p0, fmpz_poly_t pk_p1, fmpz_poly_t rlk_r0, fmpz_poly_t rlk_r1);
