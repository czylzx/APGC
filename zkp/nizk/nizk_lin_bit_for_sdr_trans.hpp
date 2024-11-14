
#ifndef NIZK_LIN_BIT_HPP_
#define NIZK_LIN_BIT_HPP_

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"
#include "../../pke/twisted_exponential_elgamal.hpp"
#include "../../commitment/pedersen.hpp"
#include "../../utility/polymul.hpp"


namespace LinBit{

using Serialization::operator<<; 
using Serialization::operator>>; 

struct PP
{
    std::vector<ECPoint> vec_g;
    ECPoint h;
    
};

// structure of instance
struct Instance
{
    ECPoint P;
};

// structure of witness 
struct Witness
{
    std::vector<BigInt> vec_a;
    BigInt r;
};


// structure of proof 
struct Proof
{
    ECPoint A;
    ECPoint C;
    ECPoint D;
    std::vector<BigInt> vec_a;
    std::vector<BigInt> vec_f;
    BigInt zA;
    BigInt zC;
    
};
 

std::ofstream &operator<<(std::ofstream &fout, const Proof &proof)
{
    fout << proof.A << proof.C << proof.D << proof.vec_a << proof.vec_f << proof.zA << proof.zC;
    return fout;  
}

std::ifstream &operator>>(std::ifstream &fin, Proof &proof)
{
    fin >> proof.A >> proof.C >> proof.D >> proof.vec_a >> proof.vec_f >> proof.zA >> proof.zC;
    return fin; 
} 

/* Setup algorithm */ 
PP Setup(size_t N)
{
    PP pp; 
    pp.vec_g= GenRandomECPointVector(N);
    pp.h = GenRandomECPoint();

    return pp; 
}
  

Proof Prove(PP &pp, Instance &instance, Witness &witness, std::string &transcript_str)
{
    size_t N = pp.vec_g.size();

    Proof proof;

    BigInt rA = GenRandomBigIntLessThan(order);
    BigInt rC = GenRandomBigIntLessThan(order);
    BigInt rD = GenRandomBigIntLessThan(order);
    std::vector<BigInt> vec_alpha = GenRandomBigIntVectorLessThan(N,order);

    //compute A
    proof.A = pp.vec_g[0] * vec_alpha[0];
    for(auto j=1; j < N; j++){
        proof.A += pp.vec_g[j] * vec_alpha[j];
    }
    proof.A = proof.A + (pp.h * rA);

    //compute C
    std::vector<BigInt> temp_C = GenRandomBigIntVectorLessThan(N,order);
    for(auto j=0; j < N; j++){
        temp_C[j] = vec_alpha[j] * (bn_1 - (witness.vec_a[j] + witness.vec_a[j]));
    }
    proof.C = pp.vec_g[0] * temp_C[0];
    for(auto j=1; j < N; j++){
        proof.C += pp.vec_g[j] * temp_C[j];
    }
    proof.C = proof.C + (pp.h * rC);

    //compute D
    std::vector<BigInt> temp_D = GenRandomBigIntVectorLessThan(N,order);
    for(auto j=0; j < N; j++){
        temp_D[j] = (-vec_alpha[j]) * vec_alpha[j];
    }
    proof.D += pp.vec_g[0] * temp_D[0];
    for(auto j=1; j < N; j++){
        proof.D += pp.vec_g[j] * temp_D[j];
    }
    proof.D = proof.D + (pp.h * rD);

    transcript_str += proof.A.ToByteString(); 
    transcript_str += proof.C.ToByteString();
    transcript_str += proof.D.ToByteString(); 

    // computer the challenge
    BigInt x = Hash::StringToBigInt(transcript_str); // V's challenge in Zq: apply FS-transform to generate the challenge

    proof.vec_f = GenRandomBigIntVectorLessThan(N,order);
    for(auto j=0; j < N; j++){
        proof.vec_f[j] = x * witness.vec_a[j] + vec_alpha[j];
    }

    proof.zA = (witness.r * x + rA) % order; 
    proof.zC = (rC * x + rD) % order;

    proof.vec_a = GenRandomBigIntVectorLessThan(N,order);
    for(auto i=0;i<N;i++){
        proof.vec_a[i] = vec_alpha[i];
    }
    
    return proof;

}


// check NIZK proof PI for Ci = Enc(pki, m; r) the witness is (r1, r2, m)
bool Verify(PP &pp, Instance &instance, std::string &transcript_str, Proof &proof)
{
    size_t N = pp.vec_g.size();
    
    transcript_str += proof.A.ToByteString();
    transcript_str += proof.C.ToByteString();
    transcript_str += proof.D.ToByteString(); 

    // compute the challenge
    BigInt x = Hash::StringToBigInt(transcript_str); // apply FS-transform to generate the challenge
    
    std::vector<bool> vec_condition(2, false);
    // check condition 1
    
    ECPoint temp1,temp2;
    temp1 = pp.vec_g[0] * proof.vec_f[0];
    for(auto j=1; j < N; j++){
        temp1 += pp.vec_g[j] * proof.vec_f[j];
    }
    temp1 = temp1 + (pp.h * proof.zA);

    temp2 = instance.P * x + proof.A;
    
    if(temp1 == temp2){vec_condition[0] = true;}

    
    // check condition 2
    ECPoint temp3,temp4;

    temp3 = pp.vec_g[0] * (proof.vec_f[0] * (x - proof.vec_f[0]));
    for(auto j=1; j < N; j++){
        temp3 += pp.vec_g[j] * (proof.vec_f[j] * (x - proof.vec_f[j]));
    }
    temp3 = temp3 + (pp.h * proof.zC);

    temp4 = (proof.C * x) + proof.D;

    if(temp3 == temp4){vec_condition[1] = true;}


    bool Validity = vec_condition[0] && vec_condition[1];


    #ifdef DEBUG
    for(auto i = 0; i < 2; i++){
        std::cout << std::boolalpha << "Condition "<< std::to_string(i) <<" (Lin Bit proof) = " 
                  << vec_condition[i] << std::endl; 
    }

    if (Validity){ 
        std::cout << "NIZK proof for Lin Bit accepts >>>" << std::endl; 
    } else {
        std::cout << "NIZK proof for Lin Bit rejects >>>" << std::endl; 
    }
    #endif

    return Validity;
}

}

#endif



