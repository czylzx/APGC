#define DEBUG

#include "../pke/twisted_exponential_elgamal.hpp"
#include "../zkp/nizk/nizk_multi_plaintext_equality.hpp"
#include "../crypto/setup.hpp"

void GenRandomInstanceWitness(MutiliPlaintextEquality::PP &pp, MutiliPlaintextEquality::Instance &instance, 
                                       MutiliPlaintextEquality::Witness &witness, bool flag)
{
    PrintSplitLine('-');  
    if (flag == true){
        std::cout << "generate n well-formed twisted ElGamal ciphertext >>>" << std::endl; 
    } else{
        std::cout << ">>> generate an ill-formed recipient twisted ElGamal ciphertext" << std::endl; 
    }

    witness.vec_r = GenRandomBigIntVectorLessThan(pp.VECTOR_LEN, order);
    witness.vec_v = GenRandomBigIntVectorLessThan(pp.VECTOR_LEN, order);
    

    instance.pk_a = GenRandomGenerator();
    instance.vec_CL.resize(pp.VECTOR_LEN);
    instance.vec_CR.resize(pp.VECTOR_LEN);
    for(auto i = 0; i < pp.VECTOR_LEN; i++)
    {
        instance.vec_CL[i] = instance.pk_a * witness.vec_r[i];
        instance.vec_CR[i] = pp.g * witness.vec_r[i] + pp.h * witness.vec_v[i];
    }
    

    if(flag == false){
        // ECPoint noisy = GenRandomGenerator();
        // instance.ct.Y = instance.ct.Y + noisy;
        for(auto i = 0; i < pp.VECTOR_LEN; i++)
        {
        instance.vec_CL[i] = GenRandomGenerator();
        instance.vec_CR[i] = GenRandomGenerator();
        }
    } 
}

void test_nizk_multi_plaintext_equality(size_t number, bool flag)
{
    PrintSplitLine('-');  
    std::cout << "begin the test of NIZKPoK for plaintext equality >>>" << std::endl; 

    TwistedExponentialElGamal::PP pp_enc = TwistedExponentialElGamal::Setup(32, 7); 
    MutiliPlaintextEquality::PP pp = MutiliPlaintextEquality::Setup(pp_enc.g, pp_enc.h, number);
    MutiliPlaintextEquality::Instance instance; 
    MutiliPlaintextEquality::Witness witness; 
    MutiliPlaintextEquality::Proof proof; 

    std::string transcript_str; 

    GenRandomInstanceWitness(pp, instance, witness, flag); 
    auto start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = ""; 
    MutiliPlaintextEquality::Prove(pp, instance, witness, transcript_str, proof); 
    auto end_time = std::chrono::steady_clock::now(); // end to count the time
    auto running_time = end_time - start_time;
    std::cout << "proof generation takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = ""; 
    MutiliPlaintextEquality::Verify(pp, instance, transcript_str, proof);
    end_time = std::chrono::steady_clock::now(); // end to count the time
    running_time = end_time - start_time;
    std::cout << "proof verification takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;
}

int main()
{
    CRYPTO_Initialize();   
    
    size_t number = 8;
    test_nizk_multi_plaintext_equality(number, true);
    test_nizk_multi_plaintext_equality(number, false); 
 
    CRYPTO_Finalize(); 

    return 0; 
}



