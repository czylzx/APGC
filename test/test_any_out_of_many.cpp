#define DEBUG

#include "../zkp/nizk/nizk_any_out_of_many.hpp"
#include "../crypto/setup.hpp"
#include "../commitment/pedersen.hpp"


void test_any_out_of_many(size_t N_max, size_t N_sender)
{
    PrintSplitLine('-');  
    std::cout << "begin the test of any_out_of_many.hpp >>>" << std::endl; 

    AnyOutOfMany::PP pp;
   
    AnyOutOfMany::Instance instance;
    AnyOutOfMany::Witness witness;
    AnyOutOfMany::Proof proof;

    Pedersen::PP pp_com;
    pp_com = Pedersen::Setup(N_max);

    pp=AnyOutOfMany::Setup(N_max, pp_com);
    
    std::vector<BigInt> vec_s(N_sender);
    std::vector<BigInt> vec_b(N_max);

    size_t sum_b = 0;
    for(auto i=0;i<N_max;i++)
    {
        srand(time(0));
        size_t random_b = rand() % 2;
        if(random_b == 0)
        {
            vec_b[i] = bn_0;
        }
        else
        {
            if(sum_b < N_sender)
            {
                vec_b[i] = bn_1;
                sum_b++;
            }
            else
            {
                vec_b[i] = bn_0;
            }      
        }  
    }

    for(auto i=0;i<N_sender;i++)
    {
        vec_s[i]=GenRandomBigIntLessThan(order);
    }
    witness.vec_s=vec_s;
    witness.vec_b=vec_b;
    std::vector<BigInt> value ={bn_0};
    size_t j = 0;
    for(auto i=0;i<N_sender;i++)
    {
        if(vec_b[i] == bn_1)
        {
            instance.vec_com.push_back(Pedersen::Commit(pp_com, value, vec_s[j]));
            j++;
        }
        else
        {
            BigInt random_value = GenRandomBigIntLessThan(order);
            instance.vec_com.push_back(Pedersen::Commit(pp_com, value, random_value));
        }
        
    }

    
   
    std::string transcript_str = "";
    AnyOutOfMany::Prove(pp, instance, witness, proof, transcript_str);
    
    transcript_str = "";
    bool testval = AnyOutOfMany::Verify(pp, instance, proof, transcript_str);
    
    if(testval == true)
    {
        std::cout<<"verify success"<<std::endl;
    }
    else
    {
        std::cout<<"verify fail"<<std::endl;
    }

}
int main()
{
    CRYPTO_Initialize(); 
    size_t N_max = 8;
    size_t N_sender;
   
    test_any_out_of_many(N_max, N_sender);

    CRYPTO_Finalize(); 

    return 0; 
}