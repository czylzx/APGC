#define DEBUG

#include "../zkp/nizk/nizk_1oon.hpp"
#include "../crypto/setup.hpp"

void GenRandomDDHInstanceWitness(_1oon::PP &pp, _1oon::Instance &instance, 
                                 _1oon::Witness &witness, bool flag)
{
    // generate a true statement (false with overwhelming probability)
    PrintSplitLine('-'); 
    if (flag == true){
        std::cout << "generate a One Out Of Many tuple >>>" << std::endl;
    } 
    else{
        std::cout << "generate a random tuple >>>" << std::endl; 
    }

    size_t m = pp.vec_g.size();
    size_t N = 2;
    for(auto i=1;i<m;i++){
        N = N * 2;
    }
     
    srand(time(0));
    witness.l = rand() % N;
    witness.r = GenRandomBigIntLessThan(order);
    
    instance.vec_c = GenRandomECPointVector(N);
    instance.vec_c[witness.l] = pp.g * witness.r; 
    
    if(flag == false){
        ECPoint noisy = GenRandomGenerator();
        instance.vec_c[witness.l] += noisy;
    }
    
}

void test_nizk_1oon(bool flag)
{
    PrintSplitLine('-');
    std::cout << "begin the test of 1oon proof >>>" << std::endl;

    size_t N_max = 8;
    
    _1oon::PP pp = _1oon::Setup(N_max);
    _1oon::Instance instance; 
    _1oon::Witness witness;  
    std::string transcript_str;


    GenRandomDDHInstanceWitness(pp, instance, witness, flag); 
    auto start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = "";
    _1oon::Proof proof = _1oon::Prove(pp, instance, witness, transcript_str);
    auto end_time = std::chrono::steady_clock::now(); // end to count the time
    auto running_time = end_time - start_time;
    std::cout << "1oon proof generation takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;


    start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = "";
    _1oon::Verify(pp, instance, transcript_str, proof);
    end_time = std::chrono::steady_clock::now(); // end to count the time
    running_time = end_time - start_time;
    std::cout << "1oon proof verification takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    std::cout << "finish the test of 1oon proof >>>" << std::endl; 

}

int main()
{
    CRYPTO_Initialize();  
    
    test_nizk_1oon(true);
    test_nizk_1oon(false); 

    CRYPTO_Finalize(); 

    return 0; 
}



