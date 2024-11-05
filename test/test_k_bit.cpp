#define DEBUG

#include "../zkp/nizk/nizk_kbit.hpp"
#include "../crypto/setup.hpp"

std::vector<size_t> Decompose(size_t l, size_t n, size_t m)
{
    std::vector<size_t> vec_index(m); 
    for(auto j = 0; j < m; j++){
        vec_index[j] = l % n;  
        l = l / n; 
    }
    return vec_index;  
} 

void GenRandomInstanceWitness(Kbit::PP &pp, Kbit::Instance &instance, 
                                 Kbit::Witness &witness, bool flag)
{
    // generate a true statement (false with overwhelming probability)
    PrintSplitLine('-'); 
    if (flag == true){
        std::cout << "generate a Kbit tuple >>>" << std::endl;
    } 
    else{
        std::cout << "generate a random tuple >>>" << std::endl; 
    }

    size_t  N = pp.vec_g.size();

    srand(time(0));
    size_t temp = rand() % N; 

    witness.r = GenRandomBigIntLessThan(order);
    std::vector<size_t> vec_coin = Decompose(temp, 2, N);
    witness.vec_a.resize(N);
    witness.vec_b.resize(N);

    size_t k = 0;

    for(auto i=0;i<N;i++){
        if(vec_coin[i] == 1){
            k++;
            witness.vec_a[i] = bn_1;
            witness.vec_b[i] = bn_0;
        }
        else{
            witness.vec_a[i] = bn_0;
            witness.vec_b[i] = bn_1.ModNegate(order);
        }
    }
    
    instance.P = pp.u * witness.r;
    for(auto i=0; i<N; i++){
        instance.P += pp.vec_g[i] * witness.vec_a[i];
        instance.P += pp.vec_h[i] * witness.vec_b[i];
    }

    instance.k = BigInt(k); 

    if(flag == false){
        ECPoint noisy = GenRandomGenerator();
        instance.P += noisy;
    } 
}

void test_nizk_LogBit(bool flag)
{
    PrintSplitLine('-');
    std::cout << "begin the test of Log Bit proof >>>" << std::endl; 
    
    size_t N_max = 8;

    Kbit::PP pp = Kbit::Setup(N_max);
    Kbit::Instance instance; 
    Kbit::Witness witness;  
    std::string transcript_str;

    GenRandomInstanceWitness(pp, instance, witness, flag); 

    auto start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = "";
    Kbit::Proof proof = Kbit::Prove(pp, instance, witness, transcript_str); 
    auto end_time = std::chrono::steady_clock::now(); // end to count the time
    auto running_time = end_time - start_time;
    std::cout << "Kbit proof generation takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;


    start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = "";
    Kbit::Verify(pp, instance, transcript_str, proof);
    end_time = std::chrono::steady_clock::now(); // end to count the time
    running_time = end_time - start_time;
    std::cout << "Kbit proof verification takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    std::cout << "finish the test of Kbit proof >>>" << std::endl; 

}

int main()
{
    CRYPTO_Initialize();  
    
    test_nizk_LogBit(true);
    test_nizk_LogBit(false); 

    CRYPTO_Finalize(); 

    return 0; 
}



