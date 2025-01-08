// #define DEBUG

#include "../zkp/nizk/nizk_ACF_1oon_range.hpp"
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

void GenRandomInstanceWitness(ACF_1oon_range::PP &pp, ACF_1oon_range::Instance &instance,
                                ACF_1oon_range::Witness &witness, bool flag, size_t k)
{

    #ifdef DEBUG
    PrintSplitLine('-'); 
    if (flag == true){
        std::cout << "generate a ACF 1oon range tuple >>>" << std::endl;
    } 
    else{
        std::cout << "generate a random tuple >>>" << std::endl; 
    }
    #endif 

    size_t n = pp.VECTOR_LEN;

    BigInt v_max = bn_1;
    for(auto i=0;i<32;i++){
        v_max = v_max * bn_2;
    }

    srand(time(0));
    witness.l = rand() % n;
    witness.vec_x = GenRandomBigIntVectorLessThan(2, order);
    witness.vec_x[1] = witness.vec_x[1] % v_max; // in [0,v_max)

    instance.vec_P = GenRandomECPointVector(n);
    instance.vec_P[witness.l] = pp.g * witness.vec_x[0] + pp.h * witness.vec_x[1];

    if(flag == false){
        ECPoint noise = GenRandomECPoint();
        instance.vec_P[witness.l] += noise;
    }
}

void test_ACF_1oon_range(bool flag)
{
    size_t VECTOR_LEN = 4;

    ACF_1oon_range::PP pp = ACF_1oon_range::Setup(VECTOR_LEN);
    ACF_1oon_range::Instance instance;
    ACF_1oon_range::Witness witness;
    std::string transcript_str = "";
    
    GenRandomInstanceWitness(pp, instance, witness,flag, 0);

    auto start_time = std::chrono::steady_clock::now(); // start to count the time

    ACF_1oon_range::Proof proof = ACF_1oon_range::Prove(pp, instance, witness, transcript_str);
    
    auto end_time = std::chrono::steady_clock::now(); // end to count the time
    auto running_time = end_time - start_time;
    std::cout << "ACF 1oon range proof generation takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;


    transcript_str = "";
    auto start_time1 = std::chrono::steady_clock::now(); // start to count the time
    
    bool result = ACF_1oon_range::Verify(pp, instance, transcript_str, proof);
    
    auto end_time1 = std::chrono::steady_clock::now(); // end to count the time
    auto running_time1 = end_time1 - start_time1;
    std::cout << "ACF 1oon range proof verify takes time = " 
    << std::chrono::duration <double, std::milli> (running_time1).count() << " ms" << std::endl;
    std::cout << result << std::endl;

}

void performance_test(bool flag)
{
    PrintSplitLine('-');
    std::cout << "begin the time test of ACF 1oon range proof >>>" << std::endl;

    std::vector<size_t> all_N = {4, 8, 16, 32, 64, 128};
    size_t round = 30;

    for(auto i=0;i<all_N.size();i++){
        size_t N = all_N[i];
        size_t K = 1;

        long long Prove_time = 0; 
        long long Verify_time = 0; 
        double proof_size = 0;

        for(auto d=0;d<round;d++){
            ACF_1oon_range::PP pp = ACF_1oon_range::Setup(N);
            ACF_1oon_range::Instance instance;
            ACF_1oon_range::Witness witness;
            std::string transcript_str = "";                
            
            GenRandomInstanceWitness(pp, instance, witness,flag, K);
                
            auto start_time = std::chrono::steady_clock::now(); 
            ACF_1oon_range::Proof proof = ACF_1oon_range::Prove(pp, instance, witness, transcript_str);
            auto end_time = std::chrono::steady_clock::now(); 
            auto running_time_prove = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            Prove_time += running_time_prove.count();

            std::ofstream fout; 
            std::ifstream fin; 
            std::string ACF_1oon_range_File = "ACF1oonRange";
            // clear file
            std::ofstream file(ACF_1oon_range_File, std::ios::trunc);
            file.close();
            // save proof
            fout.open(ACF_1oon_range_File, std::ios::binary|std::ios::app); 
            fout << proof;
            fout.close();
            // compute size
            fin.open(ACF_1oon_range_File, std::ios::ate | std::ios::binary);
            auto size = fin.tellg();
            fin.close();
            // accumulate
            proof_size += double(size);

            transcript_str = "";    
            start_time = std::chrono::steady_clock::now(); 
            bool result = ACF_1oon_range::Verify(pp, instance, transcript_str, proof);
            end_time = std::chrono::steady_clock::now();
            auto running_time_verify = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            Verify_time += running_time_verify.count();
        }
        std::cout<<"N = "<<N<<" -- "<<"K = "<<K<<std::endl;
        std::cout << "ACF 1oon range proof generation takes time = " 
        << std::chrono::duration <double, std::milli> (Prove_time).count()/round << " ms" << std::endl;
        std::cout << "ACF 1oon range proof verification takes time = " 
        << std::chrono::duration <double, std::milli> (Verify_time).count()/round << " ms" << std::endl;
        std::cout<< "ACF 1oon range proof size = "<< proof_size/round <<std::endl;
        PrintSplitLine('-');

    }

    std::cout << "finish the time test of ACF 1oon range proof >>>" << std::endl; 

}

int main()
{
    CRYPTO_Initialize();

    // test_ACF_1oon_range(true);
    // test_ACF_1oon_range(false);
    performance_test(true);

    CRYPTO_Finalize();

    return 0;
}
