#define DEBUG

#include "../zkp/nizk/nizk_log_bit.hpp"
#include "../crypto/setup.hpp"

void GenRandomInstanceWitness(LogBit::PP &pp, LogBit::Instance &instance, 
                                 LogBit::Witness &witness, bool flag)
{
    // generate a true statement (false with overwhelming probability)
    PrintSplitLine('-'); 
    if (flag == true){
        std::cout << "generate a LogBit tuple >>>" << std::endl;
    } 
    else{
        std::cout << "generate a random tuple >>>" << std::endl; 
    }

    size_t  N = pp.vec_g.size();

    witness.r = GenRandomBigIntLessThan(order);
    witness.vec_a = GenRandomBigIntVectorLessThan(N,order);
    witness.vec_b = GenRandomBigIntVectorLessThan(N,order);
    for(auto i=0;i<N;i++){
        witness.vec_b[i] = (witness.vec_a[i] - bn_1) % order;
    }

    instance.P = pp.u * witness.r;
    for(auto i=0; i<N; i++){
        instance.P += pp.vec_g[i] * witness.vec_a[i];
        instance.P += pp.vec_h[i] * witness.vec_b[i];
    }

    srand(time(0));
    instance.k = rand() % N;  

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

    LogBit::PP pp = LogBit::Setup(N_max);
    LogBit::Instance instance; 
    LogBit::Witness witness;  
    std::string transcript_str;

    GenRandomInstanceWitness(pp, instance, witness, flag); 

    auto start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = "";
    LogBit::Proof proof = LogBit::Prove(pp, instance, witness, transcript_str); 
    auto end_time = std::chrono::steady_clock::now(); // end to count the time
    auto running_time = end_time - start_time;
    std::cout << "LogBit proof generation takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;


    start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = "";
    LogBit::Verify(pp, instance, transcript_str, proof);
    end_time = std::chrono::steady_clock::now(); // end to count the time
    running_time = end_time - start_time;
    std::cout << "LogBit proof verification takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    std::cout << "finish the test of LogBit proof >>>" << std::endl; 

}

int main()
{
    CRYPTO_Initialize();  
    
    test_nizk_LogBit(true);
    // test_nizk_LogBit(false); 

    CRYPTO_Finalize(); 

    return 0; 
}



