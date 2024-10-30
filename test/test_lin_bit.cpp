#define DEBUG

#include "../zkp/nizk/nizk_lin_bit.hpp"
#include "../crypto/setup.hpp"

void GenRandomInstanceWitness(LinBit::PP &pp, LinBit::Instance &instance, 
                                 LinBit::Witness &witness, bool flag)
{
    // generate a true statement (false with overwhelming probability)
    PrintSplitLine('-'); 
    if (flag == true){
        std::cout << "generate a LinBit tuple >>>" << std::endl;
    } 
    else{
        std::cout << "generate a random tuple >>>" << std::endl; 
    }

    size_t  N = pp.vec_g.size();

    witness.r = GenRandomBigIntLessThan(order);
    witness.vec_a = GenRandomBigIntVectorLessThan(N,order);
    for(auto i=0;i<N;i++){
        witness.vec_a[i] = bn_1;
    }

    for(auto j=0; j<N; j++){
        instance.P += pp.vec_g[j] * witness.vec_a[j];
    }
    instance.P += pp.h * witness.r;

    if(flag == false){
        ECPoint noisy = GenRandomGenerator();
        instance.P = instance.P + noisy;
    } 
}

void test_nizk_LinBit(bool flag)
{
    PrintSplitLine('-');
    std::cout << "begin the test of lin bit proof >>>" << std::endl; 
    
    size_t N_max = 8;
    ECPoint h = GenRandomGenerator();
    std::vector<ECPoint> vec_g = GenRandomECPointVector(N_max);
    LinBit::PP pp = LinBit::Setup(vec_g, h, N_max);
    LinBit::Instance instance; 
    LinBit::Witness witness;  
    std::string transcript_str;

    GenRandomInstanceWitness(pp, instance, witness, flag); 
    auto start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = "";
    LinBit::Proof proof = LinBit::Prove(pp, instance, witness, transcript_str); 
    auto end_time = std::chrono::steady_clock::now(); // end to count the time
    auto running_time = end_time - start_time;
    std::cout << "LinBit proof generation takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;


    start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = "";
    LinBit::Verify(pp, instance, transcript_str, proof);
    end_time = std::chrono::steady_clock::now(); // end to count the time
    running_time = end_time - start_time;
    std::cout << "LinBit proof verification takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    std::cout << "finish the test of LinBit proof >>>" << std::endl; 

}

int main()
{
    CRYPTO_Initialize();  
    
    test_nizk_LinBit(true);
    test_nizk_LinBit(false); 

    CRYPTO_Finalize(); 

    return 0; 
}



