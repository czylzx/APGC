#define DEBUG

#include "../zkp/nizk/nizk_solvent.hpp"
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


void GenRandomInstanceWitness(Solvent::PP &pp, Solvent::Instance &instance, 
                                 Solvent::Witness &witness, bool flag)
{
    #ifdef DEBUG
    PrintSplitLine('-'); 
    if (flag == true){
        std::cout << "generate a Solvent tuple >>>" << std::endl;
    } 
    else{
        std::cout << "generate a random tuple >>>" << std::endl; 
    }
    #endif 

    size_t N = pp.N;
    size_t m = log2(N);

    // compute v_max
    BigInt v_max = bn_1;
    for(auto i=0;i<32;i++){
        v_max = v_max * bn_2;
    }

    // set witness
    srand(time(0));
    witness.l = rand() % N;
    witness.v = GenRandomBigIntLessThan(v_max);
    witness.sk = GenRandomBigIntLessThan(order);
    witness.vl = -GenRandomBigIntLessThan(v_max);
    witness.rl = GenRandomBigIntLessThan(order);

    // set instance
    instance.vec_pk = GenRandomECPointVector(N);
    instance.vec_C_L = GenRandomECPointVector(N);
    instance.vec_C_R = GenRandomECPointVector(N);
    instance.vec_CL = GenRandomECPointVector(N);
    instance.vec_CR = GenRandomECPointVector(N);
    instance.vec_C = GenRandomECPointVector(N);

    BigInt vec_random = GenRandomBigIntLessThan(order);

    instance.vec_pk[witness.l] = pp.g * witness.sk;
    instance.vec_C_R[witness.l] = pp.g * vec_random;
    instance.vec_CR[witness.l] = pp.g * witness.rl + pp.h * witness.vl;

    instance.vec_C_L[witness.l] = pp.g * ((vec_random + witness.rl) * witness.sk % order);
    instance.vec_CL[witness.l] = pp.h * ((witness.vl - witness.v) * witness.sk % order);

    if(flag == false){
        ECPoint noise = GenRandomECPoint();
        instance.vec_pk[witness.l] += noise;
    }

}

void test_nizk_Solvent(bool flag)
{
    #ifdef DEBUG
    PrintSplitLine('-'); 
    if (flag == true){
        std::cout << "generate a SdrTrans tuple >>>" << std::endl;
    } 
    else{
        std::cout << "generate a random tuple >>>" << std::endl; 
    }
    #endif 
    
    size_t N_max = 128;

    Solvent::PP pp = Solvent::Setup(N_max);
    Solvent::Instance instance; 
    Solvent::Witness witness;  
    std::string transcript_str;

    GenRandomInstanceWitness(pp, instance, witness, flag); 
    
    auto start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = "";
    Solvent::Proof proof = Solvent::Prove(pp, instance, witness, transcript_str); 
    auto end_time = std::chrono::steady_clock::now(); // end to count the time
    auto running_time = end_time - start_time;
    std::cout << "Solvent proof generation takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;


    start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = "";
    Solvent::Verify(pp, instance, transcript_str, proof);
    end_time = std::chrono::steady_clock::now(); // end to count the time
    running_time = end_time - start_time;
    std::cout << "Solvent proof verification takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    std::cout << "finish the test of Solvent proof >>>" << std::endl; 

}

void performance_test(bool flag)
{
    PrintSplitLine('-');
    std::cout << "begin the time test of Solvent proof >>>" << std::endl;

    std::vector<size_t> all_N = {4, 8, 16, 32, 64, 128};
    size_t round = 30;

    for(auto i=0;i<all_N.size();i++){
        size_t N = all_N[i];
        size_t K = 1;

        long long Prove_time = 0; 
        long long Verify_time = 0; 
        double proof_size = 0;

        for(auto d=0;d<round;d++){
            Solvent::PP pp = Solvent::Setup(N);
            Solvent::Instance instance;
            Solvent::Witness witness;
            std::string transcript_str = "";                
            
            GenRandomInstanceWitness(pp, instance, witness, flag); 

                
            auto start_time = std::chrono::steady_clock::now(); 
            Solvent::Proof proof = Solvent::Prove(pp, instance, witness, transcript_str); 
            auto end_time = std::chrono::steady_clock::now(); 
            auto running_time_prove = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            Prove_time += running_time_prove.count();

            std::ofstream fout; 
            std::ifstream fin; 
            std::string File_name = "Solvent";
            // clear file
            std::ofstream file(File_name, std::ios::trunc);
            file.close();
            // save proof
            fout.open(File_name, std::ios::binary|std::ios::app); 
            fout << proof;
            fout.close();
            // compute size
            fin.open(File_name, std::ios::ate | std::ios::binary);
            auto size = fin.tellg();
            fin.close();
            // accumulate
            proof_size += double(size);

            transcript_str = "";    
            start_time = std::chrono::steady_clock::now(); 
            Solvent::Verify(pp, instance, transcript_str, proof);
            end_time = std::chrono::steady_clock::now();
            auto running_time_verify = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            Verify_time += running_time_verify.count();
        }
        std::cout<<"N = "<<N<<" -- "<<"K = "<<K<<std::endl;
        std::cout << "Solvent proof generation takes time = " 
        << std::chrono::duration <double, std::milli> (Prove_time).count()/round << " ms" << std::endl;
        std::cout << "Solvent proof verification takes time = " 
        << std::chrono::duration <double, std::milli> (Verify_time).count()/round << " ms" << std::endl;
        std::cout<< "Solvent proof size = "<< proof_size/round <<std::endl;
        PrintSplitLine('-');

    }

    std::cout << "finish the time test of Solvent proof >>>" << std::endl; 

}

int main()
{
    CRYPTO_Initialize();  


    test_nizk_Solvent(true);
    // test_nizk_Solvent(false); 

    // performance_test(true);

    CRYPTO_Finalize(); 

    return 0; 

}
