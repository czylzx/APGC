#define DEBUG

#include "../zkp/nizk/nizk_koon.hpp"
#include "../crypto/setup.hpp"

void GenRandomInstanceWitness(Koon::PP &pp, Koon::Instance &instance, 
                                 Koon::Witness &witness, bool flag)
{
    // generate a true statement (false with overwhelming probability)
    PrintSplitLine('-'); 
    if (flag == true){
        std::cout << "generate a Koon tuple >>>" << std::endl;
    } 
    else{
        std::cout << "generate a random tuple >>>" << std::endl; 
    }

    size_t N = pp.vec_g.size();
    size_t m = log2(N);
    //pp.k = 1;

    witness.vec_r = GenRandomBigIntVectorLessThan(pp.k, order);
    witness.vec_l.resize(pp.k);


    // srand(time(0));
    //witness.vec_l[0] = rand() % N;
    witness.vec_l[0]=0;
    std::cout <<"witness.vec_l[0] =" << witness.vec_l[0] << std::endl;
    // std::cout<<pp.k<<std::endl;

    // size_t count = 1;
    // while(count < pp.k){
    //     size_t temp = rand() % N;
    //     bool flag = true;
    //     for(auto i=0;i<count;i++){
    //         if(witness.vec_l[i] == temp){
    //             flag = false;
    //             break;
    //         }
    //     }
    //     if(flag == true){
    //         witness.vec_l[count] = temp;
    //         count++;
    //     }
    // }

    // for(auto i=0;i<pp.k;i++){
    //     std::cout<<witness.vec_l[i]<<std::endl;
    // }
    
    instance.vec_c = GenRandomECPointVector(N);
    for(auto i=0;i<pp.k;i++){
        instance.vec_c[witness.vec_l[i]] = pp.g * witness.vec_r[i];
    }
    
    if(flag == false){
        ECPoint noisy = GenRandomGenerator();
        instance.vec_c[witness.vec_l[0]] += noisy;
    }
    
}

void test_nizk_koon(bool flag)
{
    PrintSplitLine('-');
    std::cout << "begin the test of Koon proof >>>" << std::endl;

    size_t N_max = 8;
    
    Koon::PP pp = Koon::Setup(N_max);
    Koon::Instance instance; 
    Koon::Witness witness;  
    std::string transcript_str;


    GenRandomInstanceWitness(pp, instance, witness, flag);
    
    auto start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = "";
    Koon::Proof proof = Koon::Prove(pp, instance, witness, transcript_str);
    auto end_time = std::chrono::steady_clock::now(); // end to count the time
    auto running_time = end_time - start_time;
    std::cout << "Koon proof generation takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = "";
    Koon::Verify(pp, instance, transcript_str, proof);
    end_time = std::chrono::steady_clock::now(); // end to count the time
    running_time = end_time - start_time;
    std::cout << "Koon proof verification takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    std::cout << "finish the test of Koon proof >>>" << std::endl; 

}

int main()
{
    CRYPTO_Initialize();  
    
    test_nizk_koon(true);
    // test_nizk_koon(false); 

    CRYPTO_Finalize(); 

    return 0; 
}



