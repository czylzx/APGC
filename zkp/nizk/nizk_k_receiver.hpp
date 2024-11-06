
#ifndef NIZK_K_RECEIVER_HPP_
#define NIZK_K_RECEIVER_HPP_

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"
#include "../../pke/twisted_exponential_elgamal.hpp"
#include "../../commitment/pedersen.hpp"
#include "../../utility/polymul.hpp"
#include "../../zkp/nizk/nizk_lin_bit.hpp"
#include "../../zkp/nizk/nizk_koon_for_kreceiver.hpp"
#include "../../zkp/bulletproofs/bullet_proof.hpp"

namespace Kreceiver{

using Serialization::operator<<; 
using Serialization::operator>>; 


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


struct PP
{
    ECPoint g;
    ECPoint h;
    ECPoint u;
    size_t k;

    // for koon
    std::vector<ECPoint> vec_g_koon;
    std::vector<ECPoint> vec_h_koon;
    std::vector<ECPoint> vec_g_mk_koon;
    ECPoint u_new_koon; // for log bit

    // for bullet
    std::vector<ECPoint> vec_g_range1;
    std::vector<ECPoint> vec_g_range2;
    std::vector<ECPoint> vec_h_range1;
    std::vector<ECPoint> vec_h_range2;
    
};

// structure of instance
struct Instance
{
    std::vector<ECPoint> vec_C;
};

// structure of witness 
struct Witness
{
    std::vector<size_t> vec_L;
    std::vector<BigInt> vec_V;
    std::vector<BigInt> vec_R;

};


// structure of proof 
struct Proof
{
    Koon::Proof koon_proof;
    Bullet::Proof bullet_proof_one;
    Bullet::Proof bullet_proof_two;

};
 

/* Setup algorithm */ 
PP Setup(size_t N)
{
    PP pp; 

    pp.h = GenRandomGenerator();
    pp.g = GenRandomGenerator();
    pp.u = GenRandomGenerator();

    srand(time(0));
    pp.k = rand() % N;
    // pp.k = 2;

    pp.u_new_koon = GenRandomGenerator();
    pp.vec_g_koon = GenRandomECPointVector(N);
    pp.vec_h_koon = GenRandomECPointVector(N);
    pp.vec_g_mk_koon = GenRandomECPointVector(log2(N)*pp.k);

    pp.vec_g_range1 = GenRandomECPointVector(32*pp.k);
    pp.vec_g_range2 = GenRandomECPointVector(32*pp.k);
    pp.vec_h_range1 = GenRandomECPointVector(32*pp.k);
    pp.vec_h_range2 = GenRandomECPointVector(32*pp.k);

    return pp; 
}


Proof Prove(PP &pp, Instance &instance, Witness &witness, std::string &transcript_str)
{    
    Proof proof;

    size_t N = pp.vec_g_koon.size();
    size_t M = log2(N);
    size_t K = pp.k;
    size_t MK = M*K;

    // koon proof
    // koon pp
    Koon::PP koon_pp = Koon::Setup(N);
    koon_pp.g = pp.g;
    koon_pp.h = pp.h;
    koon_pp.k = pp.k;
    koon_pp.u = pp.u;
    koon_pp.u_new = pp.u_new_koon;
    koon_pp.vec_g = pp.vec_g_koon;
    koon_pp.vec_h = pp.vec_h_koon;
    koon_pp.vec_g_mk = pp.vec_g_mk_koon;

    // koon instance
    Koon::Instance koon_instance;
    koon_instance.vec_c = instance.vec_C;

    // koon witness
    Koon::Witness koon_witness;
    koon_witness.vec_l = witness.vec_L;
    koon_witness.vec_r = witness.vec_R;

    // koon transcript str
    std::string koon_transcript_str = "";

    // call koon
    proof.koon_proof = Koon::Prove(koon_pp,koon_instance,koon_witness,koon_transcript_str);

    // compute e in koon
    std::string str0 = "";
    str0 += proof.koon_proof.A.ToByteString();
    str0 += proof.koon_proof.B.ToByteString();
    str0 += proof.koon_proof.C.ToByteString();
    str0 += proof.koon_proof.D.ToByteString();
    str0 += proof.koon_proof.P.ToByteString();

    BigInt e = Hash::StringToBigInt(str0);

    // compute x0 in koon
    std::string str1 = "";
    str1 += proof.koon_proof.A.ToByteString();
    str1 += proof.koon_proof.C.ToByteString();
    str1 += proof.koon_proof.D.ToByteString();
    
    BigInt x0 = Hash::StringToBigInt(str1);

    // compute vec_x^m
    std::vector<BigInt> exp_x(M+1);
    exp_x[0] = bn_1;
    for(auto i=1;i<=M;i++){
        exp_x[i] = exp_x[i-1] * x0 % order;
    }
    
    // compute tau
    BigInt tau = proof.koon_proof.zG * exp_x[M].ModInverse(order) % order;

    // range proof
    size_t range_len = 32;
    size_t max_agg_num = K;
    
    // bullet proof for v
    // bullet one pp 
    Bullet::PP bullet_pp_one = Bullet::Setup(range_len,max_agg_num);
    bullet_pp_one.g = pp.g;
    bullet_pp_one.h = pp.h;
    bullet_pp_one.u = pp.u;
    bullet_pp_one.vec_g = pp.vec_g_range1;
    bullet_pp_one.vec_h = pp.vec_h_range1;

    // bullet one instance
    Bullet::Instance bullet_instance_one;
    bullet_instance_one.C = GenRandomECPointVector(1);

    // bullet one witness
    Bullet::Witness bullet_witness_one;
    bullet_witness_one.r.resize(1);
    bullet_witness_one.r[0] = tau;
    bullet_witness_one.v = witness.vec_V;

    // bullet one transcript str
    std::string bulllet_transcript_str_one = "";

    // call bullet proof
    Bullet::Prove(bullet_pp_one,bullet_instance_one,bullet_witness_one,bulllet_transcript_str_one,proof.bullet_proof_one);
    std::cout<<"here"<<std::endl;

    // bullet proof for v-1
    // bullet two pp 
    Bullet::PP bullet_pp_two = Bullet::Setup(range_len,max_agg_num);
    bullet_pp_two.g = pp.g;
    bullet_pp_two.h = pp.h;
    bullet_pp_two.u = pp.u;
    bullet_pp_two.vec_g = pp.vec_g_range2;
    bullet_pp_two.vec_h = pp.vec_h_range2;

    // bullet one instance
    Bullet::Instance bullet_instance_two;
    bullet_instance_two.C = GenRandomECPointVector(1);

    // bullet one witness
    Bullet::Witness bullet_witness_two;
    bullet_witness_two.r.resize(1);
    bullet_witness_two.r[0] = tau;
    bullet_witness_two.v.resize(K);
    for(auto i=0;i<K;i++){
        bullet_witness_two.v[i] = (witness.vec_V[i] - bn_1 + order) % order;
    }

    // bullet one transcript str
    std::string bulllet_transcript_str_two = "";

    // call bullet proof
    Bullet::Prove(bullet_pp_two,bullet_instance_two,bullet_witness_two,bulllet_transcript_str_two,proof.bullet_proof_two);

    return proof;

}

bool Verify(PP &pp, Instance &instance, std::string &transcript_str, Proof &proof)
{

    size_t N = pp.vec_g_koon.size();
    size_t M = log2(N);
    size_t K = pp.k;
    size_t MK = M*K;

    // compute e in koon
    std::string str0 = "";
    str0 += proof.koon_proof.A.ToByteString();
    str0 += proof.koon_proof.B.ToByteString();
    str0 += proof.koon_proof.C.ToByteString();
    str0 += proof.koon_proof.D.ToByteString();
    str0 += proof.koon_proof.P.ToByteString();

    BigInt e = Hash::StringToBigInt(str0);

    // compute x0 in koon
    std::string str1 = "";
    str1 += proof.koon_proof.A.ToByteString();
    str1 += proof.koon_proof.C.ToByteString();
    str1 += proof.koon_proof.D.ToByteString();
    
    BigInt x0 = Hash::StringToBigInt(str1);

    // compute vec_x^m
    std::vector<BigInt> exp_x(M+1);
    exp_x[0] = bn_1;
    for(auto i=1;i<=M;i++){
        exp_x[i] = exp_x[i-1] * x0 % order;
    }
    
    // compute p_i,j(x)
    std::vector<std::vector<BigInt>> P_ij;
    for(auto i=0;i<K;i++){
        //compute p_j(x)
        std::vector<BigInt> vec_P = GenRandomBigIntVectorLessThan(N,order); 

        for(auto j = 0; j < N; j++){
            vec_P[j] = bn_1;
            std::vector<size_t> vec_index = Decompose(j, 2, M); 
    
            for(auto b = 0; b < M; b++){        
                if(vec_index[b] == 1){
                    vec_P[j] = vec_P[j] * proof.koon_proof.vec_f[i*M+b] % order;
                }
                else{
                    vec_P[j] = vec_P[j] * ((x0 - proof.koon_proof.vec_f[i*M+b] + order) % order) % order;
                }    
            } 
        }
        P_ij.emplace_back(vec_P);
    }

    // compute vec_e_k
    std::vector<BigInt> vec_e_k = GenBigIntPowerVector(K, e);

    // compute G
    ECPoint G_l ;
    G_l.SetInfinity();
    for(auto j=0;j<N;j++){
        BigInt index = bn_0;
        for(auto i=0;i<K;i++){
            index = (index + vec_e_k[i] * P_ij[i][j] % order) % order;
        }
        G_l += instance.vec_C[j] * index;
    }

    ECPoint G_r = proof.koon_proof.vec_C_c[0] * (- exp_x[0]);
    for(auto d=1;d<M;d++){
        G_r += proof.koon_proof.vec_C_c[d] * (- exp_x[d]);
    }

    ECPoint G = G_l + G_r;

    std::vector<bool> vec_condition(3, false);

    // check condition 1
    // koon verify
    // koon pp
    Koon::PP koon_pp = Koon::Setup(N);
    koon_pp.g = pp.g;
    koon_pp.h = pp.h;
    koon_pp.k = pp.k;
    koon_pp.u = pp.u;
    koon_pp.u_new = pp.u_new_koon;
    koon_pp.vec_g = pp.vec_g_koon;
    koon_pp.vec_h = pp.vec_h_koon;
    koon_pp.vec_g_mk = pp.vec_g_mk_koon;

    // koon instance
    Koon::Instance koon_instance;
    koon_instance.vec_c = instance.vec_C;

    // koon transcript str
    std::string koon_transcript_str = "";

    // call koon
    vec_condition[0] = Koon::Verify(koon_pp,koon_instance,koon_transcript_str,proof.koon_proof);
    

    // range proof
    size_t range_len = 32;
    size_t max_agg_num = K;

    // check condition 2    
    // bullet verify for v
    // bullet one pp 
    Bullet::PP bullet_pp_one = Bullet::Setup(range_len,max_agg_num);
    bullet_pp_one.g = pp.g;
    bullet_pp_one.h = pp.h;
    bullet_pp_one.u = pp.u;
    bullet_pp_one.vec_g = pp.vec_g_range1;
    bullet_pp_one.vec_h = pp.vec_h_range1;

    // bullet one instance
    Bullet::Instance bullet_instance_one;
    bullet_instance_one.C = GenRandomECPointVector(1);
    bullet_instance_one.C[0] = G * exp_x[M].ModInverse(order);

    // bullet one transcript str
    std::string bulllet_transcript_str_one = "";

    // call bullet proof
    vec_condition[1] = Bullet::Verify(bullet_pp_one,bullet_instance_one,bulllet_transcript_str_one,proof.bullet_proof_one);



    // check condition 3
    // bullet proof for v-1
    // bullet two pp 
    Bullet::PP bullet_pp_two = Bullet::Setup(range_len,max_agg_num);
    bullet_pp_two.g = pp.g;
    bullet_pp_two.h = pp.h;
    bullet_pp_two.u = pp.u;
    bullet_pp_two.vec_g = pp.vec_g_range1;
    bullet_pp_two.vec_h = pp.vec_h_range1;

    // bullet one instance
    ECPoint temp = pp.h * vec_e_k[0];
    for(auto i=1;i<K;i++){
        temp += pp.h * vec_e_k[i];
    }
    Bullet::Instance bullet_instance_two;
    bullet_instance_two.C = GenRandomECPointVector(1);
    bullet_instance_two.C[0] = (G * exp_x[M].ModInverse(order)) + temp.Invert();

    // bullet one transcript str
    std::string bulllet_transcript_str_two = "";

    // call bullet proof
    vec_condition[2] = Bullet::Verify(bullet_pp_two,bullet_instance_two,bulllet_transcript_str_two,proof.bullet_proof_two);


    bool Validity = vec_condition[0] && vec_condition[1] && vec_condition[2];


    #ifdef DEBUG
    for(auto i = 0; i < 3; i++){
        std::cout << std::boolalpha << "Condition "<< std::to_string(i) <<" (Kreceiver proof) = " 
                  << vec_condition[i] << std::endl; 
    }

    if (Validity){ 
        std::cout << "NIZK proof for Kreceiver accepts >>>" << std::endl; 
    } else {
        std::cout << "NIZK proof for Kreceiver rejects >>>" << std::endl; 
    }
    #endif

    return Validity;
}

}

#endif



