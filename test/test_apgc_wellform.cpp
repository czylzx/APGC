// #define DEBUG

#include "../pke/twisted_exponential_elgamal.hpp"
#include "../zkp/apgcproofs/apgc_wellform.hpp"
#include "../crypto/setup.hpp"

void GenRandomWellformInstanceWitness(WellFormProduct::PP &pp, WellFormProduct::Instance &instance,
                                      WellFormProduct::Witness &witness)
{
    PrintSplitLine('-');
    std::cout << "generate a valid twisted elgamal ciphertext >>>" << std::endl;

    witness.vec_r = GenRandomBigIntVectorLessThan(pp.VECTOR_LEN, BigInt(order));
    witness.vec_v = GenRandomBigIntVectorLessThan(pp.VECTOR_LEN, BigInt(order));

    instance.vec_pk = GenRandomECPointVector(pp.VECTOR_LEN);
    // instance.vec_CL = GenRandomECPointVector(pp.VECTOR_LEN);
    instance.vec_CL = ECPointVectorProduct(instance.vec_pk, witness.vec_r);

    std::vector<ECPoint> vec_gr(pp.VECTOR_LEN);
    for (auto i = 0; i < pp.VECTOR_LEN; i++)
    {
        vec_gr[i] = pp.g * witness.vec_r[i];
    }
    std::vector<ECPoint> vec_hv(pp.VECTOR_LEN);
    for (auto i = 0; i < pp.VECTOR_LEN; i++)
    {
        vec_hv[i] = pp.h * witness.vec_v[i];
    }
    instance.vec_CR = ECPointVectorAdd(vec_gr, vec_hv);
}

void test_apgc_wellform()
{
    // std::cout << "begin the test of NIZKPoK for plaintext knowledge >>>" << std::endl;

    TwistedExponentialElGamal::PP pp_enc = TwistedExponentialElGamal::Setup(32, 7);
    size_t VECTOR_LEN = 2;
    //WellFormProduct::PP pp = WellFormProduct::Setup(pp_enc.g, pp_enc.h, VECTOR_LEN);
    WellFormProduct::PP pp = WellFormProduct::Setup(pp_enc.g, pp_enc.h, VECTOR_LEN);
    WellFormProduct::Instance instance;
    WellFormProduct::Witness witness;

    GenRandomWellformInstanceWitness(pp, instance, witness);

    // WellFormProduct::Proof proof;
    std::string transcript_str = "";

    WellFormProduct::Proof proof;
    WellFormProduct::Prove(pp, instance, witness, transcript_str, proof);
    // auto end_time = std::chrono::steady_clock::now(); // end to count the time
    // auto running_time = end_time - start_time;
    // std::cout << "proof generation takes time = "
    // << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    // start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = "";

    bool result = WellFormProduct::Verify(pp, instance, transcript_str, proof);
    // end_time = std::chrono::steady_clock::now(); // end to count the time
    // running_time = end_time - start_time;
    std::cout << result << std::endl;
    // << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;
}

int main()
{
    CRYPTO_Initialize();

    test_apgc_wellform();

    CRYPTO_Finalize();

    return 0;
}
