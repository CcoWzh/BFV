#include "flint/fmpz_vec.h"
#include "flint/fmpz_poly.h"

extern double PI;
extern int bit_num;
extern int sigma;       // 离散高斯分布标准差
extern int B;           // [-B, B] 近似离散高斯分布
extern int d;           // 同态加密f = x^d + 1
extern int lamda;       // 系统安全参数
extern int Max;         // 以概率ro 输出x的方法：随机生成 【0， 10000】 之间的数z，如果z小于（小于等于）ro*10000则输出x，否则不输出
extern int delta;       // 系统安全参数
extern int l;
extern fmpz_t q;
extern fmpz_t neg_Delta;
extern fmpz_t t;
extern fmpz_t Delta_2;
extern fmpz_t Delta;
extern fmpz_t T;

// 设置参数,可达到128bit安全
extern void System_Param();
// 抽样
long SampleZ(int seed);
// 生成（近似）离散高斯分布
void SampleD(fmpz_poly_t v, int seed);
// R2上的多项式(均匀分布)
void Gen_R2(fmpz_poly_t v, int seed);
void Gen_Rq_div_2(fmpz_t t, int seed);
// Rq上的多项式(均匀分布)
void Gen_Rq(fmpz_poly_t v);