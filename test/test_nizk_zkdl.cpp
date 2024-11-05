#define DEBUG

#include "../zkp/nizk/nizk_zkdl.hpp"
#include "../crypto/setup.hpp"

// generate a random instance-witness pair
void GenRandomEqmdlProductInstanceWitness(ZkdlProduct::PP &pp, ZkdlProduct::Instance &instance, ZkdlProduct::Witness &witness)
{

    std::cout << "generate random (instance, witness) pair >>>" << std::endl;

    // witness.vec_a = {bn_1,bn_1,bn_1,bn_1};
    witness.vec_a = GenRandomBigIntVectorLessThan(pp.VECTOR_LEN, BigInt(order));
    instance.G = ECPointVectorMul(pp.vec_g, witness.vec_a);
}

void test_eqmdlproduct_proof()
{
    // PrintSplitLine('-');
    // std::cout << "begin the test of eqmdlproduct proof >>>" << std::endl;

    size_t VECTOR_LEN = 4;

    ZkdlProduct::PP pp = ZkdlProduct::Setup(VECTOR_LEN, true);

    ZkdlProduct::Instance instance;
    ZkdlProduct::Witness witness;

    GenRandomEqmdlProductInstanceWitness(pp, instance, witness);

    ZkdlProduct::Proof proof;

    // auto start_time = std::chrono::steady_clock::now(); // start to count the time
    std::string transcript_str = "";
    // transcript_str += instance.P.ToByteString();

    ZkdlProduct::Prove(pp, instance, witness, transcript_str, proof);

    // auto end_time = std::chrono::steady_clock::now(); // end to count the time
    // auto running_time = end_time - start_time;
    // std::cout << "eqmdlproduct proof generation takes time = "
    // << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    // start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = "";
    // transcript_str += instance.P.ToByteString();
    // pp = ZkdlProduct::Setup(VECTOR_LEN, true);
    bool va = ZkdlProduct::Verify(pp, instance, transcript_str, proof);
    // end_time = std::chrono::steady_clock::now(); // end to count the time

    std::cout << va << std::endl;

    // running_time = end_time - start_time;
    // std::cout << "eqmdlproduct proof verification takes time = "
    // << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    // start_time = std::chrono::steady_clock::now(); // start to count the time
    // transcript_str = "";
    // transcript_str += instance.P.ToByteString();
    // ZkdlProduct::FastVerify(pp, instance, transcript_str, proof);
    // end_time = std::chrono::steady_clock::now(); // end to count the time
    // running_time = end_time - start_time;
    // std::cout << "fast eqmdlproduct proof verification takes time = "
    // << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    std::cout << "finish the test of zkdlproduct proof >>>" << std::endl;
}

int main()
{
    CRYPTO_Initialize();

    test_eqmdlproduct_proof();

    CRYPTO_Finalize();

    return 0;
}