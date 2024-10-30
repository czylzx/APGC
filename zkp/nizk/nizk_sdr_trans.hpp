
#ifndef NIZK_SDR_TRANS_HPP_
#define NIZK_SDR_TRANS_HPP_

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"
#include "../../pke/twisted_exponential_elgamal.hpp"
#include "../../commitment/pedersen.hpp"
#include "../../utility/polymul.hpp"
#include "../../zkp/nizk/nizk_1oon_for_sdr_trans.hpp"
#include "../../zkp/bulletproofs/bullet_proof_for_sdr_trans.hpp"


namespace SdrTrans{

using Serialization::operator<<; 
using Serialization::operator>>; 

struct PP
{
    std::vector<ECPoint> vec_g_1oon;
    std::vector<ECPoint> vec_g_range1;
    std::vector<ECPoint> vec_g_range2;
    std::vector<ECPoint> vec_h1;
    std::vector<ECPoint> vec_h2;
    ECPoint g,h,u;
    
};

// structure of instance
struct Instance
{
    std::vector<ECPoint> vec_C;
    ECPoint B;
};

// structure of witness 
struct Witness
{
    size_t l;
    BigInt v;
    BigInt rL;
    BigInt rR;
};


// structure of proof 
struct Proof
{
    ECPoint A,C,D;
    std::vector<ECPoint> vec_C;
    ECPoint V1,V2,S1,S2;
    ECPoint T11,T12,T21,T22;
    std::vector<BigInt> vec_f;
    BigInt zA,zC,zG;
    BigInt mu1,mu2,t1,t2,tau1,tau2;
    std::vector<ECPoint> vec_L1, vec_L2, vec_R1, vec_R2;
    BigInt a1,a2,b1,b2; 
};
 

/* Setup algorithm */ 
PP Setup(size_t N)
{
    PP pp; 
    pp.vec_g_1oon = GenRandomECPointVector(log2(N));
    pp.vec_g_range1 = GenRandomECPointVector(32);
    pp.vec_g_range2 = GenRandomECPointVector(32);
    pp.vec_h1 = GenRandomECPointVector(32);
    pp.vec_h2 = GenRandomECPointVector(32);
    pp.g = generator; 
    pp.h = Hash::StringToECPoint(pp.g.ToByteString()); 
    pp.u = GenRandomECPoint();

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

    size_t m = pp.vec_g_1oon.size();
    size_t N = 2;
    for(auto i=1;i<m;i++){
        N = N * 2;
    }


    std::vector<size_t> l_index = Decompose(witness.l,2,m);

    _1oon::PP pp_1oon = _1oon::Setup(N);
    _1oon::Instance instance_1oon;
    _1oon::Witness witness_1oon;
    std::string transcript_1oon;
    //pp
    pp_1oon.g = pp.g;
    pp_1oon.u = pp.u;
    pp_1oon.h = pp.h;
    for(auto i=0;i<m;i++){
        pp_1oon.vec_g[i] = pp.vec_g_1oon[i];
    }
    
    //instance
    instance_1oon.B = instance.B;
    instance_1oon.vec_c = GenRandomECPointVector(N);
    for(auto i=0;i<N;i++){
        instance_1oon.vec_c[i] = instance.vec_C[i];
    }
    //witness
    witness_1oon.r = witness.rL;
    witness_1oon.l = witness.l;
    witness_1oon.rB = witness.rR;
    
    _1oon::Proof LB_proof = _1oon::Prove(pp_1oon, instance_1oon, witness_1oon, transcript_1oon);
    

    // 1oon proof
    // step 1
    proof.A = LB_proof.A;
    proof.C = LB_proof.C;
    proof.D = LB_proof.D;
    proof.vec_C = GenRandomECPointVector(m);
    for(auto i=0;i<m;i++){
        proof.vec_C[i] = LB_proof.vec_C[i];
    }
    // step 3
    proof.vec_f = GenRandomBigIntVectorLessThan(m,order);
    for(auto i=0;i<m;i++){
        proof.vec_f[i] = LB_proof.vec_f[i];
    }
    proof.zA = LB_proof.zA;
    proof.zC = LB_proof.zC;
    BigInt zG = LB_proof.zd;
    proof.zG = zG;
    


    std::string _1oon_str = "";
    _1oon_str += proof.A.ToByteString();
    _1oon_str += proof.C.ToByteString();
    _1oon_str += proof.D.ToByteString();
    // for(auto i=0;i<m;i++){
    //     _1oon_str += proof.vec_C[i].ToByteString();
    // }

    // computer the challenge x_0
    BigInt x0 = Hash::StringToBigInt(_1oon_str);

    // compute x_0 ^ m
    BigInt x0m = bn_1;
    for(auto i=0;i<m;i++){
        x0m = x0m * x0 % order;
    }
    
    // compute tau
    BigInt tau = zG * x0m.ModInverse(order) % order;
    
    // range proof
    size_t range_len = 32;
    size_t max_agg_num = 1;
    
    // Bullet proof v
    Bullet::PP First_Bullet_pp = Bullet::Setup(range_len,max_agg_num);
    Bullet::Instance First_Bullet_instance;
    Bullet::Witness First_Bullet_witness;
    std::string First_Bullet_str;
    Bullet::Proof First_Bullet_proof;
    // set pp
    First_Bullet_pp.g = pp.g;
    First_Bullet_pp.h = pp.h;
    First_Bullet_pp.u = pp.u;
    First_Bullet_pp.vec_g = pp.vec_g_range1;
    First_Bullet_pp.vec_h = pp.vec_h1;
    
    // set instance
    First_Bullet_instance.C = GenRandomECPointVector(max_agg_num);
    // set witness
    First_Bullet_witness.r = GenRandomBigIntVectorLessThan(max_agg_num,order);
    First_Bullet_witness.v = GenRandomBigIntVectorLessThan(max_agg_num,order);
    First_Bullet_witness.r[0] = tau;
    First_Bullet_witness.v[0] = witness.v;
    // First_Bullet_witness.r[0] = witness.v;
    // First_Bullet_witness.v[0] = tau;

    // set transcript_str
    First_Bullet_str = "";
    
    // call first Bullet proof

    Bullet::Prove(First_Bullet_pp, First_Bullet_instance, First_Bullet_witness, First_Bullet_str, First_Bullet_proof);


    // second range proof
    Bullet::PP Second_Bullet_pp = Bullet::Setup(range_len,max_agg_num);
    Bullet::Instance Second_Bullet_instance;
    Bullet::Witness Second_Bullet_witness;
    std::string Second_Bullet_str = "";
    Bullet::Proof Second_Bullet_proof;

    // set Bullet proof v-1
    // set pp
    Second_Bullet_pp.g = pp.g;
    Second_Bullet_pp.h = pp.h;
    Second_Bullet_pp.u = pp.u;
    Second_Bullet_pp.vec_g = pp.vec_g_range2;
    Second_Bullet_pp.vec_h = pp.vec_h2;
    // set instance
    Second_Bullet_instance.C = GenRandomECPointVector(max_agg_num);
    // set witness
    Second_Bullet_witness.r = GenRandomBigIntVectorLessThan(max_agg_num,order);
    Second_Bullet_witness.v = GenRandomBigIntVectorLessThan(max_agg_num,order);
    Second_Bullet_witness.r[0] = tau;
    Second_Bullet_witness.v[0] = witness.v - bn_1;

    // set transcript_str
    Second_Bullet_str = "";

    // call Second bullet proof
    Bullet::Prove(Second_Bullet_pp, Second_Bullet_instance, Second_Bullet_witness, Second_Bullet_str, Second_Bullet_proof);

    // step 1
    proof.V1 = First_Bullet_proof.A;
    proof.S1 = First_Bullet_proof.S;
    proof.V2 = Second_Bullet_proof.A;
    proof.S2 = Second_Bullet_proof.S;

    // step 2
    proof.T11 = First_Bullet_proof.T1;
    proof.T12 = First_Bullet_proof.T2;
    proof.T21 = Second_Bullet_proof.T1;
    proof.T22 = Second_Bullet_proof.T2;

    // set transcript str
    std::string first_str = "";
    std::string second_str = "";

    // compute y1,z1
    first_str += proof.V1.ToByteString();
    BigInt y1 = Hash::StringToBigInt(first_str);

    first_str += proof.S1.ToByteString();
    BigInt z1 = Hash::StringToBigInt(first_str);

    // compute y2,z2
    second_str += proof.V2.ToByteString();
    BigInt y2 = Hash::StringToBigInt(second_str);

    second_str += proof.S2.ToByteString();
    BigInt z2 = Hash::StringToBigInt(second_str);
    
    // compute x1,x2
    first_str += proof.T11.ToByteString() + proof.T12.ToByteString();
    second_str += proof.T21.ToByteString() + proof.T22.ToByteString();
    BigInt x1 = Hash::StringToBigInt(first_str);
    BigInt x2 = Hash::StringToBigInt(second_str);

    // step 3
    proof.mu1 = First_Bullet_proof.mu;
    proof.t1 = First_Bullet_proof.tx;
    proof.vec_L1 = First_Bullet_proof.ip_proof.vec_L;
    proof.vec_R1 = First_Bullet_proof.ip_proof.vec_R;
    proof.mu2 = Second_Bullet_proof.mu;
    proof.t2 = Second_Bullet_proof.tx;
    proof.vec_L2 = Second_Bullet_proof.ip_proof.vec_L;
    proof.vec_R2 = Second_Bullet_proof.ip_proof.vec_R;
    

    // set tau1,tau2
    proof.tau1 = First_Bullet_proof.taux;
    proof.tau2 = Second_Bullet_proof.taux;

    // set a,b
    proof.a1 = First_Bullet_proof.ip_proof.a;
    proof.b1 = First_Bullet_proof.ip_proof.b;
    proof.a2 = Second_Bullet_proof.ip_proof.a;
    proof.b2 = Second_Bullet_proof.ip_proof.b;

    std::vector<ECPoint> aaa = {pp.g * zG,pp.g * zG + pp.h * witness.v};
    PrintECPointVector(aaa,"");
    return proof;

}


// check NIZK proof PI for Ci = Enc(pki, m; r) the witness is (r1, r2, m)
bool Verify(PP &pp, Instance &instance, std::string &transcript_str, Proof &proof)
{
   
    size_t m = pp.vec_g_1oon.size();
    size_t N = 2;
    for(auto i=1;i<m;i++){
        N = N * 2;
    }

    std::string str = "";
    str += proof.A.ToByteString();
    str += proof.C.ToByteString();
    str += proof.D.ToByteString();
    // for(auto i=0;i<m;i++){
    //     str += proof.vec_C[i].ToByteString();
    // }

    // computer the challenge x_0
    BigInt x0 = Hash::StringToBigInt(str);

    // compute G
    std::vector<BigInt> exp_x(m+1);
    exp_x[0] = bn_1;  
    for(auto k = 1; k <= m; k++){
        exp_x[k] = exp_x[k-1] * x0 % order; 
    }

    BigInt x0m_inv = exp_x[m].ModInverse(order);  

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
                vec_P[j] = vec_P[j] * ((x0 - proof.vec_f[b] + order) % order) % order;
            }    
        } 
    }
    
    ECPoint left = instance.vec_C[0] * vec_P[0];
    for(auto i=1;i<N;i++){
        left += instance.vec_C[i] * vec_P[i];
    }

    ECPoint G = left + right;
    
    // check
    std::vector<bool> vec_condition(3, false);
    // check condition 1 -- 1oon
    _1oon::PP pp_1oon = _1oon::Setup(N);
    _1oon::Instance instance_1oon;
    _1oon::Proof proof_1oon;
    std::string transcript_1oon = "";
    // pp
    pp_1oon.g = pp.g;
    pp_1oon.u = pp.u;
    pp_1oon.h = pp.h;
    pp_1oon.vec_g = pp.vec_g_1oon;
    // instance
    instance_1oon.B = instance.B;
    instance_1oon.vec_c = instance.vec_C;
    // proof
    proof_1oon.A = proof.A;
    proof_1oon.B = instance.B;
    proof_1oon.C = proof.C;
    proof_1oon.D = proof.D;
    proof_1oon.vec_C = proof.vec_C;
    proof_1oon.vec_f = proof.vec_f;
    proof_1oon.zA = proof.zA;
    proof_1oon.zC = proof.zC;
    proof_1oon.zd = GenRandomBigIntLessThan(order);
    vec_condition[0] = _1oon::Verify(pp_1oon,instance_1oon,transcript_str,proof_1oon);


    size_t range_len = 32;
    size_t max_agg_num = 1;

    // check condition 2 -- range v
    Bullet::PP First_Bullet_pp = Bullet::Setup(range_len,max_agg_num);
    Bullet::Instance First_Bullet_instance;
    std::string First_Bullet_str;
    Bullet::Proof First_Bullet_proof;

    // set pp
    First_Bullet_pp.g = pp.g;
    First_Bullet_pp.h = pp.h;
    First_Bullet_pp.u = pp.u;
    First_Bullet_pp.vec_g = pp.vec_g_range1;
    First_Bullet_pp.vec_h = pp.vec_h1;
    // set instance
    First_Bullet_instance.C = GenRandomECPointVector(1);
    First_Bullet_instance.C[0] = G * (bn_0 - x0m_inv + order); 

    std::vector<ECPoint> aaa = {G,G};
    PrintECPointVector(aaa,"");

    // set transcript_str
    First_Bullet_str = "";
    // set proof
    First_Bullet_proof.A = proof.V1;
    First_Bullet_proof.mu = proof.mu1;
    First_Bullet_proof.S = proof.S1;
    First_Bullet_proof.T1 = proof.T11;
    First_Bullet_proof.T2 = proof.T12;
    First_Bullet_proof.taux = proof.tau1;
    First_Bullet_proof.tx = proof.t1;
    First_Bullet_proof.ip_proof.vec_L = proof.vec_L1;
    First_Bullet_proof.ip_proof.vec_R = proof.vec_R1;
    First_Bullet_proof.ip_proof.a = proof.a1;
    First_Bullet_proof.ip_proof.b = proof.b1;

    vec_condition[1] = Bullet::Verify(First_Bullet_pp,First_Bullet_instance,First_Bullet_str,First_Bullet_proof);

    // check condition 3 -- range (v-1)
    Bullet::PP Second_Bullet_pp = Bullet::Setup(range_len,max_agg_num);
    Bullet::Instance Second_Bullet_instance;
    std::string Second_Bullet_str;
    Bullet::Proof Second_Bullet_proof;

    // set pp
    Second_Bullet_pp.g = pp.g;
    Second_Bullet_pp.h = pp.h;
    Second_Bullet_pp.u = pp.u;
    Second_Bullet_pp.vec_g = pp.vec_g_range2;
    Second_Bullet_pp.vec_h = pp.vec_h2;
    // set instance
    Second_Bullet_instance.C = GenRandomECPointVector(1);
    Second_Bullet_instance.C[0] = G * (bn_0 - x0m_inv + order) + pp.h.Invert(); 


    // std::vector<BigInt> aaa = {x0}
    // set transcript_str
    Second_Bullet_str = "";
    // set proof
    Second_Bullet_proof.A = proof.V2;
    Second_Bullet_proof.mu = proof.mu2;
    Second_Bullet_proof.S = proof.S2;
    Second_Bullet_proof.T1 = proof.T21;
    Second_Bullet_proof.T2 = proof.T22;
    Second_Bullet_proof.taux = proof.tau2;
    Second_Bullet_proof.tx = proof.t2;
    Second_Bullet_proof.ip_proof.vec_L = proof.vec_L2;
    Second_Bullet_proof.ip_proof.vec_R = proof.vec_R2;
    Second_Bullet_proof.ip_proof.a = proof.a2;
    Second_Bullet_proof.ip_proof.b = proof.b2;


    vec_condition[2] = Bullet::Verify(Second_Bullet_pp,Second_Bullet_instance,Second_Bullet_str,Second_Bullet_proof);



    bool Validity = vec_condition[0] && vec_condition[1] && vec_condition[2];


    #ifdef DEBUG
    for(auto i = 0; i < 3; i++){
        std::cout << std::boolalpha << "Condition "<< std::to_string(i) <<" (Sdr Trans proof) = " 
                  << vec_condition[i] << std::endl; 
    }

    if (Validity){ 
        std::cout << "NIZK proof for Sdr Trans accepts >>>" << std::endl; 
    } else {
        std::cout << "NIZK proof for Sdr Trans rejects >>>" << std::endl; 
    }
    #endif

    return Validity;
}

}

#endif



