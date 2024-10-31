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
    // std::cout << "pp.vec_g.size = " << pp.vec_g.size() << std::endl;
    // std::cout << "vec_bit.size = " << vec_bit.size() << std::endl;
    instance.B = ECPointVectorMul(pp.vec_g, vec_bit) + pp.h * witness.rb_l0;
    instance.pk = GenRandomECPointVector(pp.VECTOR_LEN);
    instance.pk[sender_index] = pp.g * sk;
   

    size_t n = pp.VECTOR_LEN;
    std::vector<TwistedExponentialElGamal::CT> vec_participant_transfer_ct(n);
    std::vector<TwistedExponentialElGamal::CT> vec_participant_balance_ct(n);

    BigInt v_balance = BigInt(128);
    for(auto i = 0; i < n; i++)
    {
        vec_participant_balance_ct[i] = TwistedExponentialElGamal::Enc(pp_enc, instance.pk[i], v_balance, GenRandomBigIntLessThan(order));
        vec_participant_transfer_ct[i] = TwistedExponentialElGamal::Enc(pp_enc, instance.pk[i], bn_0, GenRandomBigIntLessThan(order));
    }
    BigInt v = BigInt(32);
    vec_participant_transfer_ct[sender_index] = TwistedExponentialElGamal::Enc(pp_enc, instance.pk[sender_index], -v, GenRandomBigIntLessThan(order));
    vec_participant_transfer_ct[receiver_index] = TwistedExponentialElGamal::Enc(pp_enc, instance.pk[receiver_index], v, GenRandomBigIntLessThan(order));
    vec_participant_balance_ct[sender_index] = TwistedExponentialElGamal::HomoAdd(vec_participant_balance_ct[sender_index], vec_participant_transfer_ct[sender_index]);
    vec_participant_balance_ct[receiver_index] = TwistedExponentialElGamal::HomoAdd(vec_participant_balance_ct[receiver_index], vec_participant_transfer_ct[receiver_index]);

    witness.r_refresh = GenRandomBigIntLessThan(order);
    std::vector<ECPoint> sum_ct_left(n);
    std::vector<ECPoint> sum_ct_right(n);
    for(auto i = 0; i < n; i++)
    {
        sum_ct_left[i] = vec_participant_transfer_ct[i].X + vec_participant_balance_ct[i].X;
        sum_ct_right[i] = vec_participant_transfer_ct[i].Y + vec_participant_balance_ct[i].Y;
    }
    instance.Sum_CL = sum_ct_left;
    instance.Sum_CR = sum_ct_right;
    witness.balance_sender = BigInt(96);
    instance.Refresh_CL = instance.pk[sender_index] * witness.r_refresh;
    instance.Refresh_CR = pp_enc.g * witness.r_refresh + pp_enc.h * witness.balance_sender ;
    
}

void test_apgc_solvent_equal()
{
    // std::cout << "begin the test of NIZKPoK for plaintext knowledge >>>" << std::endl;

    TwistedExponentialElGamal::PP pp_enc = TwistedExponentialElGamal::Setup(32, 7);
    size_t VECTOR_LEN = 4;

    Solvent_Equal::PP pp = Solvent_Equal::Setup(pp_enc.g, pp_enc.h, VECTOR_LEN);
    
    Solvent_Equal::Instance instance;
    Solvent_Equal::Witness witness;
    
    GenRandomInstanceWitness(pp_enc, pp, instance, witness);


    std::string transcript_str = "";

    Solvent_Equal::Proof proof;
    Solvent_Equal::Prove(pp, instance, witness, transcript_str, proof);

    transcript_str = "";

    bool result = Solvent_Equal::Verify(pp, instance, transcript_str, proof);
 
    std::cout << result << std::endl;

}

int main()
{
    CRYPTO_Initialize();

    test_apgc_solvent_equal();

    CRYPTO_Finalize();

    return 0;
}
