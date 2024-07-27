
#include "unit.h"
#include "flint/fmpz_vec.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpz.h"

void poly_to_num(fmpz_poly_t c,fmpz_t result)
{
    slong degree;  // 多项式的度数
    ulong coeff;   // 多项式的系数
    // 计算多项式对应的整数
    degree = fmpz_poly_degree(c);
    fmpz_set_ui(result, 0); // 初始化结果为0
    for (slong i = 0; i <= degree; i++)
    {
        coeff = fmpz_poly_get_coeff_ui(c, i); // 获取第i个系数
        if (coeff == 1)                       // 如果系数为1，则在结果中添加相应的二进制位
        {
            fmpz_add_ui(result, result, 1 << i); // 将结果左移i位，然后加1
        }
    }
    // 打印结果
    // printf("多项式对应的整数是: %s\n", fmpz_get_str(NULL, 10, result));
}