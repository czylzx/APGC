
#ifndef NIZK_1OON_HPP_
#define NIZK_1OON_HPP_

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"
#include "../../pke/twisted_exponential_elgamal.hpp"
#include "../../commitment/pedersen.hpp"
#include "../../utility/polymul.hpp"
#include "../../zkp/nizk/nizk_lin_bit.hpp"


namespace _1oon{

using Serialization::operator<<; 
using Serialization::operator>>; 

struct PP
{
    std::vector<ECPoint> vec_g;
    ECPoint g;
    ECPoint h;
    ECPoint u;
    
};

// structure of instance
struct Instance
{
    std::vector<ECPoint> vec_c;
};

// structure of witness 
struct Witness
{
    size_t l;
    BigInt r;
};


// structure of proof 
struct Proof
{
    ECPoint A;
    ECPoint B;
    ECPoint C;
    ECPoint D;
    std::vector<ECPoint> vec_C;
    std::vector<BigInt> vec_f;
    BigInt zA;
    BigInt zC;
    BigInt zd;
};
 

/* Setup algorithm */ 
PP Setup(size_t N)
{
    PP pp; 
    pp.vec_g = GenRandomECPointVector(logb(N));

    pp.h = GenRandomGenerator();
    pp.g = GenRandomGenerator();
    pp.u = GenRandomGenerator();

    return pp; 
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


Proof Prove(PP &pp, Instance &instance, Witness &witness, std::string &transcript_str)
{    
    Proof proof;

    size_t m = pp.vec_g.size();
    size_t N = 2;
    for(auto i=1;i<m;i++){
        N = N * 2;
    }

    std::vector<size_t> vec_l = Decompose(witness.l, 2, m);

    BigInt rB = GenRandomBigIntLessThan(order);
    
    proof.B = pp.vec_g[0] * vec_l[0];
    for(auto i=1; i<m; i++){
        proof.B += pp.vec_g[i] * vec_l[i] ;
    }
    proof.B += pp.u * rB;
    
    //proof lin bit
    LinBit::PP linbit_pp = LinBit::Setup(pp.vec_g, pp.u, m);
    /*linbit_pp.h = pp.u;
    for(auto i=0;i<m;i++){
        linbit_pp.vec_g[i] = pp.vec_g[i];
    }*/
    
    LinBit::Instance linbit_instance; 
    linbit_instance.P = proof.B;
    
    LinBit::Witness linbit_witness;
    linbit_witness.vec_a = GenRandomBigIntVectorLessThan(m,order);
    for(auto i=0; i<m; i++){
        linbit_witness.vec_a[i] = vec_l[i]; 
    }
    linbit_witness.r = rB;

    std::string transcript_str_linbit = "";
    
    LinBit::Proof LB_proof = LinBit::Prove(linbit_pp, linbit_instance, linbit_witness, transcript_str_linbit);

    proof.A = LB_proof.A;
    proof.C = LB_proof.C;
    proof.D = LB_proof.D;
    proof.zA = LB_proof.zA;
    proof.zC = LB_proof.zC;
    proof.vec_f = GenRandomBigIntVectorLessThan(m,order);// only define is enough
    for(auto i=0;i<m;i++){
        proof.vec_f[i] = LB_proof.vec_f[i];
    }
    
    std::vector<BigInt> vec_a = GenRandomBigIntVectorLessThan(m,order);
    for(auto i=0;i<m;i++){
        vec_a[i] = LB_proof.vec_a[i];
    }

    std::string str = "";
    str += proof.A.ToByteString();
    str += proof.C.ToByteString();
    str += proof.D.ToByteString();
    // computer the challenge
    BigInt x = Hash::StringToBigInt(str); // V's challenge in Zq: apply FS-transform to generate the challenge

    std::vector<BigInt> vec_rho = GenRandomBigIntVectorLessThan(m,order);
    
    //compute p_j_d
    std::vector<std::vector<BigInt>> P; 

    for(auto j = 0; j < N; j++){
        std::vector<std::vector<BigInt>> A(m, std::vector<BigInt>(2));        
        // prepare m ploynomial of form ax+b
        std::vector<size_t> vec_index = Decompose(j, 2, m); 
 
        for(auto b = 0; b < m; b++){      
            if(vec_index[b] == 1){
                A[b][0] = vec_a[b];
                A[b][1] = vec_l[b];
            }
            else{
                A[b][0] = bn_0 - vec_a[b];
                A[b][1] = bn_1 - vec_l[b];
            }    
        } 
        std::vector<BigInt> p_j = PolyMul(A);
            
        P.emplace_back(p_j); 
    }

    //compute C_deg

    proof.vec_C = GenRandomECPointVector(m);
    for(auto d=0; d < m; d++){
        proof.vec_C[d] = instance.vec_c[0] * P[0][d];
        for(auto j=1; j<N; j++){
            proof.vec_C[d] += instance.vec_c[j] * P[j][d];
        }
        proof.vec_C[d] += pp.g * vec_rho[d];
    }

    //compute zd
    std::vector<BigInt> exp_x(m+1);
    exp_x[0] = bn_1;  
    for(auto k = 1; k <= m; k++){
        exp_x[k] = exp_x[k-1] * x % order; 
    }

    proof.zd = witness.r * exp_x[m] % order;
    for(auto d=0;d<m;d++){
        proof.zd = ((proof.zd - (vec_rho[d] * exp_x[d] % order))+ order) % order;
    }

    return proof;

}

bool Verify(PP &pp, Instance &instance, std::string &transcript_str, Proof &proof)
{
    size_t m = pp.vec_g.size();
    size_t N = 2;
    for(auto i=1;i<m;i++){
        N = N * 2;
    }

    transcript_str += proof.A.ToByteString();
    transcript_str += proof.C.ToByteString();
    transcript_str += proof.D.ToByteString();
     
    // compute the challenge
    BigInt x = Hash::StringToBigInt(transcript_str); 

    std::vector<bool> vec_condition(2, false);

    // check condition 1
    LinBit::PP linbit_pp = LinBit::Setup(pp.vec_g, pp.u, m);
    /*linbit_pp.h = pp.u;
    for(auto i=0;i<m;i++){
        linbit_pp.vec_g[i] = pp.vec_g[i];
    }*/
    
    LinBit::Instance linbit_instance; 
    linbit_instance.P = proof.B;
    
    LinBit::Proof linbit_proof;
    linbit_proof.A = proof.A;
    linbit_proof.C = proof.C;
    linbit_proof.D = proof.D;
    linbit_proof.zA = proof.zA;
    linbit_proof.zC = proof.zC;
    linbit_proof.vec_f = GenRandomBigIntVectorLessThan(m,order);
    for(auto i=0;i<m;i++){
        linbit_proof.vec_f[i] = proof.vec_f[i];
    }
    
    std::string transcript_str_linbit = "";
    vec_condition[0] = LinBit::Verify(linbit_pp, linbit_instance, transcript_str_linbit, linbit_proof);
    
    // check condition 2

    std::vector<BigInt> exp_x(m+1);
    exp_x[0] = bn_1;  
    for(auto k = 1; k <= m; k++){
        exp_x[k] = exp_x[k-1] * x % order; 
    }

    //right part
    ECPoint right = proof.vec_C[0] * (bn_0 - exp_x[0] + order);
    for(auto d=1;d<m;d++){
        right += proof.vec_C[d] * (bn_0 - exp_x[d] + order);
    }
    
    //left part
    //compute p_j(x)
    std::vector<BigInt> vec_P = GenRandomBigIntVectorLessThan(N,order); 

    for(auto j = 0; j < N; j++){
        vec_P[j] = bn_1;
        std::vector<size_t> vec_index = Decompose(j, 2, m); 
 
        for(auto b = 0; b < m; b++){        
            if(vec_index[b] == 1){
                vec_P[j] = vec_P[j] * proof.vec_f[b] % order;
            }
            else{
                vec_P[j] = vec_P[j] * ((x - proof.vec_f[b] + order) % order) % order;
            }    
        } 
    }
    
    ECPoint left = instance.vec_c[0] * vec_P[0];
    for(auto i=1;i<N;i++){
        left += instance.vec_c[i] * vec_P[i];
    }

    if((left + right) == (pp.g * proof.zd)){ 
        vec_condition[1] = true; 
    }

    bool Validity = vec_condition[0] && vec_condition[1];


    #ifdef DEBUG
    for(auto i = 0; i < 2; i++){
        std::cout << std::boolalpha << "Condition "<< std::to_string(i) <<" (1oon proof) = " 
                  << vec_condition[i] << std::endl; 
    }

    if (Validity){ 
        std::cout << "NIZK proof for One out of many accepts >>>" << std::endl; 
    } else {
        std::cout << "NIZK proof for One out of many rejects >>>" << std::endl; 
    }
    #endif

    return Validity;
}

}

#endif



