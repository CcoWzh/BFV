#pragma once
#include "parameter.h"
#include "flint/fmpz_vec.h"
#include "flint/fmpz_poly.h"

// 将明文表示为二进制
void Encode(int m, fmpz_poly_t pm);
// 加密和解密
void SH_Encrypt(fmpz_poly_t m, fmpz_poly_t pk_p0, fmpz_poly_t pk_p1, fmpz_poly_t c0, fmpz_poly_t c1);
// void SH_Decrypt(fmpz_poly_t sk, fmpz_poly_t c0, fmpz_poly_t c1);
void SH_Decrypt(fmpz_poly_t sk, fmpz_poly_t c0, fmpz_poly_t c1,fmpz_poly_t dec_m);

// 执行共同的计算[\lfloor\frac{t(*)}{q}\rceil]_q
void fmpz_poly_nearest_fmpz(fmpz_poly_t v);

void SH_Add(fmpz_poly_t c_00, fmpz_poly_t c_01, fmpz_poly_t c_10, fmpz_poly_t c_11, 
            fmpz_poly_t c0, fmpz_poly_t c1);
void SH_Mul(fmpz_poly_t c_00, fmpz_poly_t c_01, fmpz_poly_t c_10, fmpz_poly_t c_11, 
            fmpz_poly_t c0, fmpz_poly_t c1, fmpz_poly_t c2);
// 对乘法后的结果进行解密，获得解密后的数据            
void SH_DecMul(fmpz_poly_t sk, fmpz_poly_t c0, fmpz_poly_t c1, fmpz_poly_t c2);
// 重线性化
void Relinear(fmpz_poly_t *rlk_r0, fmpz_poly_t *rlk_r1, fmpz_poly_t c0, fmpz_poly_t c1, fmpz_poly_t c2, 
            fmpz_poly_t _c0, fmpz_poly_t _c1);
