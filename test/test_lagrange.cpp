// #define DEBUG

#include "../pke/twisted_exponential_elgamal.hpp"
#include "../zkp/apgcproofs/apgc_amorhom_proof.hpp"
#include "../crypto/setup.hpp"
#include <vector>
#include <iostream>

std::vector<BigInt> lagrange(std::vector<BigInt> x, std::vector<BigInt> y)
    {
        size_t n = x.size();
        if(n != y.size())
        {
            std::cerr << "vector size does not match!" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::vector<BigInt> result(n, bn_0);
        for(auto i = 0; i < n; i++)
        {
            BigInt num = bn_1;
            BigInt den = bn_1;
            for(auto j = 0; j < n; j++)
            {
                if(i != j)
                {
                    num = num * x[j] % order;
                    den = den * (x[j] - x[i]) % order;
                }
            }
            result[i] = y[i] * num * den.ModInverse(order)%order;
        }
        return result;
    }

void test()
{
    size_t n = 4;
    std::vector<BigInt> vex_x(n);
    std::vector<BigInt> vec_y(n, bn_0);
    for(auto i = 0; i < n; i++)
    {
        vex_x[i] = BigInt(i);
    }
    vec_y[0] = bn_1;
    std::vector<BigInt> vec_P = lagrange(vex_x, vec_y);
    for(auto i = 0; i < n; i++)
    {
        vec_P[i].Print("P[i]");
    }
}

int main()
{
    CRYPTO_Initialize();

    test();

    CRYPTO_Finalize();

    return 0;
}
