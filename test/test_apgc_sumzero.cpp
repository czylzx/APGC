#define DEBUG

#include "../zkp/apgcproofs/apgc_sum_zero.hpp"
#include "../crypto/setup.hpp"
#include "../commitment/pedersen.hpp"

void GenRandomSumzeroProductInstanceWitness(SumZero::PP &pp, SumZero::Instance &instance, SumZero::Witness &witness)
{

    std::cout << "generate random (instance, witness) pair >>>" << std::endl;

    // witness.vec_a = {bn_1,bn_1,bn_1,bn_1};
    witness.r = GenRandomBigIntVectorLessThan(pp.n, BigInt(order));
    // instance.C = GenRandomECPoint();
    BigInt sumr = bn_0;
    for (auto i = 0; i < pp.n; i++)
    {
        sumr += witness.r[i];
    }
    instance.C = pp.g * sumr;
}

void test_sumzero_proof()
{
    PrintSplitLine('-');
    std::cout << "begin the test of sum_zero.hpp >>>" << std::endl;

    SumZero::PP pp = SumZero::Setup(4);
    SumZero::Instance instance;
    SumZero::Witness witness;
    GenRandomSumzeroProductInstanceWitness(pp, instance, witness);

    SumZero::Proof proof;
    std::string transcript_str = "";
    SumZero::Prove(pp, instance, witness, transcript_str, proof);

    bool va = SumZero::Verify(pp, instance, transcript_str, proof);

    if (va == true)
    {
        std::cout << "verify success" << std::endl;
    }
    else
    {
        std::cout << "verify fail" << std::endl;
    }
}

int main()
{
    CRYPTO_Initialize();

    test_sumzero_proof();

    CRYPTO_Finalize();

    return 0;
}