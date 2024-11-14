// #define DEBUG

#include "../zkp/compare/nizk_ACF_koon.hpp"
#include "../crypto/setup.hpp"

size_t GetTheNthBit(size_t index, size_t n)
{
    return (index>>n)&1;
}

/* generate a^n = (a^0, a^1, a^2, ..., a^{n-1}) */ 
std::vector<BigInt> GenBigIntPowerVector(size_t LEN, const BigInt &a)
{
    std::vector<BigInt> vec_result(LEN);
    vec_result[0] = BigInt(bn_1); 
    for (auto i = 1; i < LEN; i++)
    {
        vec_result[i] = (vec_result[i-1] * a) % order; // result[i] = result[i-1]*a % order
    }
    return vec_result; 
}

std::vector<size_t> Decompose(size_t l, size_t n, size_t m)
{
    std::vector<size_t> vec_index(m); 
    for(auto j = 0; j < m; j++){
        vec_index[j] = l % n;  
        l = l / n; 
    }
    return vec_index;  
} 

void GenRandomInstanceWitness(ACF_koon::PP &pp, ACF_koon::Instance &instance,
                                ACF_koon::Witness &witness, bool flag)
{
    PrintSplitLine('-'); 
    if (flag == true){
        std::cout << "generate a ACF Koon tuple >>>" << std::endl;
    } 
    else{
        std::cout << "generate a random tuple >>>" << std::endl; 
    }

    size_t n = pp.VECTOR_LEN;

    srand(time(0));
    // instance.k = rand() % n;
    instance.k = 32;
    instance.vec_c = GenRandomECPointVector(n);

    witness.vec_r = GenRandomBigIntVectorLessThan(n,order);
    witness.vec_s = GenRandomBigIntVectorLessThan(n,order);

    for(auto i=0;i<instance.k;i++){
        witness.vec_s[i] = bn_1;
    }
    for(auto i=instance.k;i<n;i++){
        witness.vec_s[i] = bn_0;
    }

    for(auto j=0;j<n;j++){
        instance.vec_c[j] = pp.g * (witness.vec_r[j] * witness.vec_s[j] % order);
    }

    if(flag == false){
        ECPoint noise = GenRandomECPoint();
        instance.vec_c[0] += noise;
    }
}

void test_ACF_koon(bool flag)
{
    size_t VECTOR_LEN = 64;

    ACF_koon::PP pp = ACF_koon::Setup(VECTOR_LEN);
    ACF_koon::Instance instance;
    ACF_koon::Witness witness;
    std::string transcript_str = "";
    
    
    GenRandomInstanceWitness(pp, instance, witness,flag);

    auto start_time = std::chrono::steady_clock::now(); // start to count the time

    ACF_koon::Proof proof = ACF_koon::Prove(pp, instance, witness, transcript_str);
    
    auto end_time = std::chrono::steady_clock::now(); // end to count the time
    auto running_time = end_time - start_time;
    std::cout << "ACF Koon proof generation takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;


    transcript_str = "";
    auto start_time1 = std::chrono::steady_clock::now(); // start to count the time
    
    bool result = ACF_koon::Verify(pp, instance, transcript_str, proof);
    
    auto end_time1 = std::chrono::steady_clock::now(); // end to count the time
    auto running_time1 = end_time1 - start_time1;
    std::cout << "ACF Koon proof verify takes time = " 
    << std::chrono::duration <double, std::milli> (running_time1).count() << " ms" << std::endl;
    std::cout << result << std::endl;

}

int main()
{
    CRYPTO_Initialize();

    test_ACF_koon(true);
    // test_ACF_koon(false);

    CRYPTO_Finalize();

    return 0;
}
