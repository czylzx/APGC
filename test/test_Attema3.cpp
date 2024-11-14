// #define DEBUG

#include "../zkp/nizk/nizk_Attema3.hpp"
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

void GenRandomInstanceWitness(Attema3::PP &pp, Attema3::Instance &instance,
                                Attema3::Witness &witness, bool flag)
{
    PrintSplitLine('-'); 
    if (flag == true){
        std::cout << "generate a Attema3 tuple >>>" << std::endl;
    } 
    else{
        std::cout << "generate a random tuple >>>" << std::endl; 
    }

    size_t n = pp.n;
    size_t k = pp.k;

    witness.vec_S.resize(k);
    
    srand(time(0));
    witness.vec_S[0] = rand() % n;
    
    size_t count = 1;
    while(count < k){
        size_t temp = rand() % n;
        bool flag = true;
        for(auto i=0;i<count;i++){
            if(witness.vec_S[i] == temp){
                flag = false;
                break;
            }
        }
        if(flag == true){
            witness.vec_S[count] = temp;
            count++;
        }
    }

    for(auto i=0;i<k;i++){
        witness.vec_S[i] += 1;
    }  

    witness.vec_x = GenRandomBigIntVectorLessThan(n,order);

    instance.vec_P.resize(n);
    for(auto i=0;i<k;i++){
        instance.vec_P[witness.vec_S[i]-1] = pp.g * witness.vec_x[witness.vec_S[i]-1];
    }

    if(flag == false){
        ECPoint noise = GenRandomECPoint();
        instance.vec_P[witness.vec_S[0]] += noise;
    }
}

void test_Attema3(size_t N, size_t K, bool flag)
{

    Attema3::PP pp = Attema3::Setup(N, K);
    Attema3::Instance instance;
    Attema3::Witness witness;
    std::string transcript_str = "";
    
    
    GenRandomInstanceWitness(pp, instance, witness, flag);

    auto start_time = std::chrono::steady_clock::now(); // start to count the time

    Attema3::Proof proof = Attema3::Prove(pp, instance, witness, transcript_str);
    
    auto end_time = std::chrono::steady_clock::now(); // end to count the time
    auto running_time = end_time - start_time;
    std::cout << "Attema3 proof generation takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;


    transcript_str = "";
    auto start_time1 = std::chrono::steady_clock::now(); // start to count the time
    
    bool result = Attema3::Verify(pp, instance, transcript_str, proof);
    
    auto end_time1 = std::chrono::steady_clock::now(); // end to count the time
    auto running_time1 = end_time1 - start_time1;
    std::cout << "Attema3 proof verify takes time = " 
    << std::chrono::duration <double, std::milli> (running_time1).count() << " ms" << std::endl;
    std::cout << result << std::endl;

}

int main()
{
    CRYPTO_Initialize();

    size_t N = 8;
    size_t K = 4;
    test_Attema3(N, K, true);
    test_Attema3(N,K,false);

    CRYPTO_Finalize();

    return 0;
}
