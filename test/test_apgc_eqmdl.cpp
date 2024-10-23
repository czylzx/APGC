#define DEBUG

#include "../zkp/apgcproofs/apgc_eqmdl_proof.hpp"
#include "../crypto/setup.hpp"

// generate a random instance-witness pair
void GenRandomEqmdlProductInstanceWitness(EqmdlProduct::PP &pp, EqmdlProduct::Instance &instance, EqmdlProduct::Witness &witness)
{ 

    std::cout << "generate random (instance, witness) pair >>>" << std::endl;  

    // InnerProduct_Instance_new(instance); 
    witness.vec_a = GenRandomBigIntVectorLessThan(pp.VECTOR_LEN, BigInt(order)); 
    instance.G = GenRandomECPoint();
    instance.H = GenRandomECPoint();
    instance.P = GenRandomECPoint();
    // // witness.vec_b = GenRandomBigIntVectorLessThan(pp.VECTOR_LEN, BigInt(order)); 

    // //instance.u = GenRandomGenerator();
    // BigInt c = BigIntVectorModInnerProduct(witness.vec_a BigInt(order)); 

    // instance.P = pp.u * c;  // P = u^c
 
    // instance.P = instance.P + ECPointVectorMul(pp.vec_g, witness.vec_a) ;
}

void test_eqmdlproduct_proof()
{
    PrintSplitLine('-');
    std::cout << "begin the test of eqmdlproduct proof >>>" << std::endl; 
    
    size_t VECTOR_LEN = 32; 

    EqmdlProduct::PP pp = EqmdlProduct::Setup(VECTOR_LEN, true);
    
    EqmdlProduct::Instance instance; 
    EqmdlProduct::Witness witness; 

    GenRandomEqmdlProductInstanceWitness(pp, instance, witness); 


    EqmdlProduct::Proof proof; 

    auto start_time = std::chrono::steady_clock::now(); // start to count the time
    std::string transcript_str = ""; 
    transcript_str += instance.P.ToByteString(); 

    EqmdlProduct::Prove(pp, instance, witness, transcript_str, proof);
    
    auto end_time = std::chrono::steady_clock::now(); // end to count the time
    auto running_time = end_time - start_time;
    std::cout << "eqmdlproduct proof generation takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = ""; 
    transcript_str += instance.P.ToByteString(); 
    pp = EqmdlProduct::Setup(VECTOR_LEN, true);
    bool va = EqmdlProduct::Verify(pp, instance, transcript_str, proof); 
    end_time = std::chrono::steady_clock::now(); // end to count the time

    std::cout << va << std::endl;

    running_time = end_time - start_time;
    std::cout << "eqmdlproduct proof verification takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;


    

    start_time = std::chrono::steady_clock::now(); // start to count the time
    transcript_str = ""; 
    transcript_str += instance.P.ToByteString(); 
    // EqmdlProduct::FastVerify(pp, instance, transcript_str, proof); 
    end_time = std::chrono::steady_clock::now(); // end to count the time
    running_time = end_time - start_time;
    std::cout << "fast eqmdlproduct proof verification takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    std::cout << "finish the test of eqmdlproduct proof >>>" << std::endl; 
}

int main()
{
    CRYPTO_Initialize();  
    
    test_eqmdlproduct_proof();

    CRYPTO_Finalize(); 

    return 0; 
}