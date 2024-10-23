//#define DEBUG

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
    instance.vec_CL = GenRandomECPointVector(pp.VECTOR_LEN); 
    instance.vec_CR = GenRandomECPointVector(pp.VECTOR_LEN); 
    // TwistedExponentialElGamal::PP pp_enc; 
    // pp_enc.g = pp.g; 
    // pp_enc.h = pp.h;  

    // instance.ct = TwistedExponentialElGamal::Enc(pp_enc, instance.pk, witness.v, witness.r); 
}

void test_apgc_wellform()
{
    std::cout << "begin the test of NIZKPoK for plaintext knowledge >>>" << std::endl; 
    
    TwistedExponentialElGamal::PP pp_enc = TwistedExponentialElGamal::Setup(32, 7); 
    size_t VECTOR_LEN = 32; 
    WellFormProduct::PP pp = WellFormProduct::Setup(VECTOR_LEN,pp_enc);
    WellFormProduct::Instance instance;
    WellFormProduct::Witness witness; 

    GenRandomWellformInstanceWitness(pp, instance, witness); 

    std::string transcript_str; 

    auto start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = ""; 
    WellFormProduct::Proof proof = WellFormProduct::Prove(pp, instance, witness, transcript_str); 
    auto end_time = std::chrono::steady_clock::now(); // end to count the time
    auto running_time = end_time - start_time;
    std::cout << "proof generation takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = ""; 
    WellFormProduct::Verify(pp, instance, transcript_str, proof); 
    end_time = std::chrono::steady_clock::now(); // end to count the time
    running_time = end_time - start_time;
    std::cout << "proof verification takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;
}

int main()
{
    CRYPTO_Initialize(); 
    
    test_apgc_wellform();

    CRYPTO_Finalize(); 

    return 0; 
}

