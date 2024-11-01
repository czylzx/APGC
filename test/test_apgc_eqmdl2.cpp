#define DEBUG

#include "../zkp/apgcproofs/apgc_eqmdl_proof.hpp"
#include "../crypto/setup.hpp"

// generate a random instance-witness pair
void GenRandomEqmdlProductInstanceWitness(EqmdlProduct::PP &pp, EqmdlProduct::Instance &instance, EqmdlProduct::Witness &witness)
{

    std::cout << "generate random (instance, witness) pair >>>" << std::endl;

    // witness.vec_a = {bn_1,bn_1,bn_1,bn_1};
    witness.vec_a = GenRandomBigIntVectorLessThan(pp.VECTOR_LEN, BigInt(order));
    instance.G = ECPointVectorMul(pp.vec_g, witness.vec_a);
    instance.P = ECPointVectorMul(pp.vec_p, witness.vec_a);
}

void test_eqmdlproduct_proof()
{
    // PrintSplitLine('-');
    // std::cout << "begin the test of eqmdlproduct proof >>>" << std::endl;

    size_t VECTOR_LEN = 4;

    EqmdlProduct::PP pp = EqmdlProduct::Setup(VECTOR_LEN, true);

    EqmdlProduct::Instance instance;
    EqmdlProduct::Witness witness;

    GenRandomEqmdlProductInstanceWitness(pp, instance, witness);

    EqmdlProduct::Proof proof;

    std::string transcript_str = "";

    EqmdlProduct::Prove(pp, instance, witness, transcript_str, proof);


    transcript_str = "";

    bool va = EqmdlProduct::Verify(pp, instance, transcript_str, proof);

    std::cout << va << std::endl;

    std::cout << "finish the test of eqmdlproduct proof >>>" << std::endl;
}

int main()
{
    CRYPTO_Initialize();

    test_eqmdlproduct_proof();

    CRYPTO_Finalize();

    return 0;
}