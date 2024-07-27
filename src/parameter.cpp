#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "parameter.h"
#include "flint/fmpz_vec.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpz.h"

double PI = 3.1415926;
int bit_num = 128;
int sigma = 8;      // 离散高斯分布标准差
int B = 10 * sigma; //[-B, B] 近似离散高斯分布
int d = 8192;       // 同态加密f = x^d + 1
int Max = 10000;    // 以概率ro 输出x的方法：随机生成 【0， 10000】 之间的数z，如果z小于（小于等于）ro*10000则输出x，否则不输出
int delta;          // 系统安全参数
int l;
// 哦，还能这么做啊!!
fmpz_t q;
fmpz_t neg_Delta;
fmpz_t t;
fmpz_t Delta_2;
fmpz_t Delta;
fmpz_t T;

void System_Param()
{
    // PI = 3.1415926;
    // Max = 10000;
    // sigma = 8;
    // B = 10 * sigma;
    // d = 8192; // 2^13 = 8192
    // bit_num = 128;
    const char *str1 = "100000000000000000000000000000000";
    const char *str2 = "8000";
    fmpz_set_str(q, str1, 16); // q = 34028236692093846346337460743176821145665536 = pow(2, 128)
    fmpz_set_str(t, str2, 16); // t = pow(2, 15)=32768
    fmpz_cdiv_q(Delta, q, t);  // Delta = |q/t|
    fmpz_neg(neg_Delta, Delta);
    fmpz_cdiv_q_ui(Delta_2, Delta, 2);
    fmpz_sqrt(T, q); // pow(2, 64) = 18446744073709551616
    l = fmpz_flog(q, T);
}

long SampleZ(int seed)
{
    long xx;
    long z;
    double ro;
    long zz;
    int b;

    srand(time(NULL) + seed);

    do
    {
        xx = rand() % B;
        z = rand() % Max;
        ro = exp(-(PI * xx * xx) / (sigma * sigma));
    } while (!(z < ro * Max));

    b = rand() % 2;
    if (b == 0)
        b = -1;
    xx = xx * b;

    return xx;
}

void SampleD(fmpz_poly_t v, int seed)
{
    int i;
    long x;

    for (i = d - 1; i >= 0; i--)
    {
        x = SampleZ(i + seed);
        fmpz_poly_set_coeff_si(v, i, x);
    }
}

void Gen_R2(fmpz_poly_t v, int seed)
{
    int i;
    srand(time(NULL) + seed);
    int b;
    for (i = d - 1; i >= 0; i--)
    {
        b = rand() % 2;
        fmpz_poly_set_coeff_ui(v, i, b);
    }
}

void Gen_Rq_div_2(fmpz_t t, int seed)
{
    int i;
    char *str = (char *)malloc(sizeof(char) * (bit_num - 1));
    int b;
    srand(time(NULL) + seed);

    for (i = bit_num - 2; i >= 0; i--)
        str[i] = rand() % 2 + 48;
    fmpz_set_str(t, str, 2);
    b = rand() % 2;

    if (b == 0)
        b = -1;
    fmpz_mul_si(t, t, b);
}

void Gen_Rq(fmpz_poly_t v)
{
    int i;
    fmpz_t t;
    fmpz_init(t);
    for (i = d - 1; i >= 0; i--)
    {
        Gen_Rq_div_2(t, i);
        fmpz_poly_set_coeff_fmpz(v, i, t);
    }

    fmpz_clear(t);
}
