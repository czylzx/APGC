#define DEBUG

#include "../zkp/nizk/nizk_sdr_trans.hpp"
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

void GenRandomInstanceWitness(SdrTrans::PP &pp, SdrTrans::Instance &instance, 
                                 SdrTrans::Witness &witness, bool flag)
{
    // generate a true statement (false with overwhelming probability)
    PrintSplitLine('-'); 
    if (flag == true){
        std::cout << "generate a SdrTrans tuple >>>" << std::endl;
    } 
    else{
        std::cout << "generate a random tuple >>>" << std::endl; 
    }

    size_t m = pp.vec_g_1oon.size();
    size_t N = 2;
    for(auto i=1;i<m;i++){
        N = N * 2;
    }

    srand(time(0));
    witness.l = rand() % N;
    witness.v = -BigInt(32);
    witness.rL = GenRandomBigIntLessThan(order);
    witness.rR = GenRandomBigIntLessThan(order);

    std::vector<size_t> l_index = Decompose(witness.l, 2, m);

    instance.B = pp.h * witness.rR;
    for(auto b=0; b<m; b++){
        instance.B +=  pp.vec_g_1oon[b] * l_index[b]; 
    }

    instance.vec_C = GenRandomECPointVector(N);
    instance.vec_C[witness.l] = pp.g * witness.rL + pp.h * witness.v;

    if(flag == false){
        ECPoint noise = GenRandomECPoint();
        instance.B += noise;
        instance.vec_C[witness.l] += noise;
    }
}

void test_nizk_SdrTrans(bool flag)
{
    PrintSplitLine('-');
    std::cout << "begin the test of Sdr Trans proof >>>" << std::endl; 
    
    size_t N_max = 8;

    SdrTrans::PP pp = SdrTrans::Setup(N_max);
    SdrTrans::Instance instance; 
    SdrTrans::Witness witness;  
    std::string transcript_str;

    GenRandomInstanceWitness(pp, instance, witness, flag); 
    auto start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = "";
    
    SdrTrans::Proof proof = SdrTrans::Prove(pp, instance, witness, transcript_str); 

    auto end_time = std::chrono::steady_clock::now(); // end to count the time
    auto running_time = end_time - start_time;
    std::cout << "SdrTrans proof generation takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;


    start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = "";
    SdrTrans::Verify(pp, instance, transcript_str, proof);
    end_time = std::chrono::steady_clock::now(); // end to count the time
    running_time = end_time - start_time;
    std::cout << "SdrTrans proof verification takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    std::cout << "finish the test of SdrTrans proof >>>" << std::endl; 

}

int main()
{
    CRYPTO_Initialize();  
    
    test_nizk_SdrTrans(true);
    //test_nizk_SdrTrans(false); 

    CRYPTO_Finalize(); 

    return 0; 
}



