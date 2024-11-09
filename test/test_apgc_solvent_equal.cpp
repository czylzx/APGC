// #define DEBUG

#include "../pke/twisted_exponential_elgamal.hpp"
#include "../zkp/apgcproofs/apgc_sdr_solvent_equal.hpp"
#include "../crypto/setup.hpp"
size_t GetTheNthBit(size_t index, size_t n)
{
    return (index>>n)&1;
}
void GenRandomInstanceWitness(TwistedExponentialElGamal::PP pp_enc, Solvent_Equal::PP &pp, Solvent_Equal::Instance &instance,
                                      Solvent_Equal::Witness &witness)
{
    PrintSplitLine('-');
    std::cout << "generate a valid twisted elgamal ciphertext >>>" << std::endl;
    size_t sender_index = 0;
    size_t receiver_index = 1;

    
    BigInt sk = GenRandomBigIntLessThan(order);
    witness.sk = sk;
    witness.l0 = sender_index;
    witness.rb_l0 = GenRandomBigIntLessThan(order);
    size_t m = log2(pp.VECTOR_LEN);

    std::vector<BigInt> vec_bit(m);
    for(auto i = 0;i < m; i++)
    {
        if(GetTheNthBit(sender_index, i) == 1)
        {
            vec_bit[i] = bn_1;
        }
        else
        {
            vec_bit[i] = bn_0;
        }
    }

    instance.B = ECPointVectorMul(pp.vec_g, vec_bit) + pp.h * witness.rb_l0;
    instance.pk = GenRandomECPointVector(pp.VECTOR_LEN);
    instance.pk[witness.l0] = pp.g * sk;
   
    witness.balance_sender = BigInt(32);
    witness.r_refresh = GenRandomBigIntLessThan(order);

    size_t n = pp.VECTOR_LEN;

    instance.Refresh_CL = instance.pk[witness.l0] * witness.r_refresh;
    instance.Refresh_CR = pp_enc.g * witness.r_refresh + pp_enc.h * witness.balance_sender ;

    instance.Sum_CL = GenRandomECPointVector(n);
    instance.Sum_CR = GenRandomECPointVector(n);

    ECPoint R = (instance.Refresh_CR - instance.Sum_CR[witness.l0]) * witness.sk;
    instance.Sum_CL[witness.l0] = instance.Refresh_CL - R;
}

void test_apgc_solvent_equal()
{
    TwistedExponentialElGamal::PP pp_enc = TwistedExponentialElGamal::Setup(32, 7);
    size_t VECTOR_LEN = 32;

    Solvent_Equal::PP pp = Solvent_Equal::Setup(pp_enc.g, pp_enc.h, VECTOR_LEN);
    Solvent_Equal::Instance instance;
    Solvent_Equal::Witness witness;
    std::string transcript_str = "";
    Solvent_Equal::Proof proof;
    
    GenRandomInstanceWitness(pp_enc, pp, instance, witness);

    auto start_time = std::chrono::steady_clock::now(); // start to count the time
    Solvent_Equal::Prove(pp, instance, witness, transcript_str, proof);
    auto end_time = std::chrono::steady_clock::now(); // end to count the time
    auto running_time = end_time - start_time;
    std::cout << "sdrsolventequal proof generation takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;


    transcript_str = "";
    auto start_time1 = std::chrono::steady_clock::now(); // start to count the time
    bool result = Solvent_Equal::Verify(pp, instance, transcript_str, proof);
    auto end_time1 = std::chrono::steady_clock::now(); // end to count the time
    auto running_time1 = end_time1 - start_time1;
    std::cout << "sdrsolventequal proof verify takes time = " 
    << std::chrono::duration <double, std::milli> (running_time1).count() << " ms" << std::endl;
    std::cout << result << std::endl;

}

int main()
{
    CRYPTO_Initialize();

    test_apgc_solvent_equal();

    CRYPTO_Finalize();

    return 0;
}
