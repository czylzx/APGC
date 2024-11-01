// #define DEBUG

#include "../pke/twisted_exponential_elgamal.hpp"
#include "../zkp/apgcproofs/apgc_amorhom_proof.hpp"
#include "../crypto/setup.hpp"

void GenRandomInstanceWitness(AmorHom::PP &pp, AmorHom::Instance &instance,
                                      AmorHom::Witness &witness)
{
    PrintSplitLine('-');
    std::cout << "generate a valid twisted elgamal ciphertext >>>" << std::endl;

    
}

void test_apgc_wellform()
{
    // std::cout << "begin the test of NIZKPoK for plaintext knowledge >>>" << std::endl;

    TwistedExponentialElGamal::PP pp_enc = TwistedExponentialElGamal::Setup(32, 7);
    size_t VECTOR_LEN = 4;

    AmorHom::PP pp = AmorHom::Setup(VECTOR_LEN);
    AmorHom::Instance instance;
    AmorHom::Witness witness;

    GenRandomInstanceWitness(pp, instance, witness);


    std::string transcript_str = "";

    AmorHom::Proof proof;
    AmorHom::Prove(pp, instance, witness, transcript_str, proof);

    transcript_str = "";

    bool result = AmorHom::Verify(pp, instance, transcript_str, proof);
 
    std::cout << result << std::endl;

}

int main()
{
    CRYPTO_Initialize();

    test_apgc_wellform();

    CRYPTO_Finalize();

    return 0;
}
