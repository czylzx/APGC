#define DEBUG

#include "../zkp/nizk/nizk_k_receiver.hpp"
#include "../crypto/setup.hpp"

void GenRandomInstanceWitness(Kreceiver::PP &pp, Kreceiver::Instance &instance, 
                                 Kreceiver::Witness &witness, bool flag)
{
    // generate a true statement (false with overwhelming probability)
    PrintSplitLine('-'); 
    if (flag == true){
        std::cout << "generate a Koon tuple >>>" << std::endl;
    } 
    else{
        std::cout << "generate a random tuple >>>" << std::endl; 
    }

    size_t N = pp.vec_g_koon.size();
    size_t m = log2(N);
    size_t k = pp.k;

    witness.vec_L.resize(k);
    witness.vec_R.resize(k);
    witness.vec_V.resize(k);

    witness.vec_L[0] = rand() % N;
    size_t count = 1;
    while(count < pp.k){
        size_t temp = rand() % N;
        bool flag = true;
        for(auto i=0;i<count;i++){
            if(witness.vec_L[i] == temp){
                flag = false;
                break;
            }
        }
        if(flag == true){
            witness.vec_L[count] = temp;
            count++;
        }
    }

    witness.vec_R = GenRandomBigIntVectorLessThan(k,order);
    witness.vec_V = GenRandomBigIntVectorLessThan(k,order);

    instance.vec_C = GenRandomECPointVector(N);
    for(auto i=0;i<k;i++){
        size_t index = witness.vec_L[i];
        instance.vec_C[index] = pp.g * witness.vec_R[i] + pp.h * witness.vec_V[i];
    }

    if(flag == false){
        ECPoint Noise = GenRandomECPoint();
        instance.vec_C[witness.vec_L[0]] += Noise;
    }
    
}

void test_nizk_Kreceiver(bool flag)
{
    PrintSplitLine('-');
    std::cout << "begin the test of Kreceiver proof >>>" << std::endl;

    size_t N_max = 8;
    
    Kreceiver::PP pp = Kreceiver::Setup(N_max);
    Kreceiver::Instance instance; 
    Kreceiver::Witness witness;  
    std::string transcript_str;

    GenRandomInstanceWitness(pp, instance, witness, flag);
    
    auto start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = "";
    Kreceiver::Proof proof = Kreceiver::Prove(pp, instance, witness, transcript_str);
    auto end_time = std::chrono::steady_clock::now(); // end to count the time
    auto running_time = end_time - start_time;
    std::cout << "Kreceiver proof generation takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = "";
    Kreceiver::Verify(pp, instance, transcript_str, proof);
    end_time = std::chrono::steady_clock::now(); // end to count the time
    running_time = end_time - start_time;
    std::cout << "Kreceiver proof verification takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    std::cout << "finish the test of Kreceiver proof >>>" << std::endl; 

}

int main()
{
    CRYPTO_Initialize();  
    
    test_nizk_Kreceiver(true);
    // test_nizk_Kreceiver(false); 

    CRYPTO_Finalize(); 

    return 0; 
}



