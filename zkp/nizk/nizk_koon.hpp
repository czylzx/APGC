
#ifndef NIZK_KOON_HPP_
#define NIZK_KOON_HPP_

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"
#include "../../pke/twisted_exponential_elgamal.hpp"
#include "../../commitment/pedersen.hpp"
#include "../../utility/polymul.hpp"
#include "../../zkp/nizk/nizk_lin_bit.hpp"
#include "../../zkp/nizk/nizk_log_bit.hpp"


namespace Koon{

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
    size_t k;
    ECPoint g;
    ECPoint h;
    ECPoint u;
    std::vector<ECPoint> vec_g;
    std::vector<ECPoint> vec_h;
    std::vector<ECPoint> vec_g_mk;
    // for log bit
    ECPoint u_new;
    
};

// structure of instance
struct Instance
{
    std::vector<ECPoint> vec_c;
};

// structure of witness 
struct Witness
{
    std::vector<size_t> vec_l;
    std::vector<BigInt> vec_r;
};


// structure of proof 
struct Proof
{
    ECPoint A;
    ECPoint B;
    ECPoint C;
    ECPoint D;
    ECPoint P;
    std::vector<ECPoint> vec_C_c;
    std::vector<ECPoint> vec_C_p;
    std::vector<BigInt> vec_f;
    BigInt zA;
    BigInt zC;
    BigInt zG;
    BigInt zP;
    LogBit::Proof logbit_proof;
};
 

/* Setup algorithm */ 
PP Setup(size_t N)
{
    PP pp; 

    srand(time(0));
    pp.k = rand() % N;
    pp.h = GenRandomGenerator();
    pp.g = GenRandomGenerator();
    pp.u = GenRandomGenerator();
    pp.u_new = GenRandomGenerator();
    pp.vec_g = GenRandomECPointVector(N);
    pp.vec_h = GenRandomECPointVector(N);
    pp.vec_g_mk = GenRandomECPointVector(log2(N) * pp.k);

    return pp; 
}

Proof Prove(PP &pp, Instance &instance, Witness &witness, std::string &transcript_str)
{    
    Proof proof;

    size_t N = pp.vec_g.size();
    size_t m = log2(N);
    size_t k = pp.k;
    size_t mk = m * k;

    // compute vec_L
    std::vector<BigInt> vec_L(mk, bn_0);

    for(auto i=0;i<k;i++){
        std::vector<size_t> vec_li = Decompose(witness.vec_l[i], 2, m);
        for(auto j=0;j<m;j++){
            vec_L[m*i+j] = vec_li[j];
        }
    }

    // compute vec_s
    std::vector<BigInt> vec_s(N, bn_0);

    for(auto i=0;i<k;i++){
        vec_s[witness.vec_l[i]] = bn_1;
    }

    // choose rB,rP
    BigInt rB = GenRandomBigIntLessThan(order);
    BigInt rP = GenRandomBigIntLessThan(order);

    // compute B
    proof.B = pp.u * rB;
    for(auto i=0;i<mk;i++){
        proof.B += pp.vec_g_mk[i] * vec_L[i];
    }

    // compute P
    proof.P = pp.u * rP;
    for(auto i=0;i<N;i++){
        proof.P += pp.vec_g[i] * vec_s[i];
        proof.P += pp.vec_h[i] * ((vec_s[i] - bn_1 + order) % order);
    }

    // Lin Bit Proof
    LinBit::PP linbit_pp = LinBit::Setup(mk);
    linbit_pp.h = pp.u;
    linbit_pp.vec_g = pp.vec_g_mk;
    
    LinBit::Instance linbit_instance; 
    linbit_instance.P = proof.B;
    
    LinBit::Witness linbit_witness;
    linbit_witness.vec_a.resize(mk);
    linbit_witness.vec_a = vec_L;
    linbit_witness.r = rB;

    std::string transcript_str_linbit = "";
    
    LinBit::Proof LB_proof = LinBit::Prove(linbit_pp, linbit_instance, linbit_witness, transcript_str_linbit);

    proof.A = LB_proof.A;
    proof.C = LB_proof.C;
    proof.D = LB_proof.D;
    proof.zA = LB_proof.zA;
    proof.zC = LB_proof.zC;
    proof.vec_f.resize(mk);
    proof.vec_f = LB_proof.vec_f;
    
    std::vector<BigInt> vec_a(mk);
    vec_a = LB_proof.vec_a;

    // computer the challenge e
    std::string str = "";
    str += proof.A.ToByteString();
    str += proof.B.ToByteString();
    str += proof.C.ToByteString();
    str += proof.D.ToByteString();
    str += proof.P.ToByteString();

    BigInt e = Hash::StringToBigInt(str);

    // choose vec_rho
    std::vector<BigInt> vec_rho = GenRandomBigIntVectorLessThan(m,order);
    
    // choose vec_alpha
    std::vector<BigInt> vec_alpha = GenRandomBigIntVectorLessThan(m,order);

    //compute p_i,j,d
    std::vector<std::vector<std::vector<BigInt>>> P;
    for(auto i=0;i<k;i++){
        std::vector<std::vector<BigInt>> Pi; 
        for(auto j = 0; j < N; j++){
            std::vector<std::vector<BigInt>> A(m, std::vector<BigInt>(2));        
            // prepare m ploynomial of form ax+b
            std::vector<size_t> vec_index = Decompose(j, 2, m); 
    
            for(auto b = 0; b < m; b++){      
                if(vec_index[b] == 1){
                    A[b][0] = vec_a[m*i+b];
                    A[b][1] = vec_L[m*i+b];
                }
                else{
                    A[b][0] = bn_0 - vec_a[m*i+b];
                    A[b][1] = bn_1 - vec_L[m*i+b];
                }    
            } 
            std::vector<BigInt> p_j = PolyMul(A);
                
            Pi.emplace_back(p_j); 
        }
        P.emplace_back(Pi);
    }

    // compute vec_e_k
    std::vector<BigInt> vec_e_k = GenBigIntPowerVector(k, e);

    //compute vec_C_c,vec_C_p
    proof.vec_C_c.resize(m);
    proof.vec_C_p.resize(m);

    for(auto d=0; d < m; d++){
        proof.vec_C_c[d] = pp.g * vec_rho[d];
        proof.vec_C_p[d] = pp.u * vec_alpha[d];

        for(auto j=0; j<N; j++){
            BigInt index1 = bn_0;
            for(auto i=0;i<k;i++){
                index1 = (index1 + vec_e_k[i] * P[i][j][d] % order) % order;
            }
            proof.vec_C_c[d] += instance.vec_c[j] * index1;

            BigInt index2 = bn_0;
            for(auto i=0;i<k;i++){
                index2 = (index2 + P[i][j][d]) % order;
            }
            proof.vec_C_p[d] += (pp.vec_g[j] + pp.vec_h[j]) * index2;
        }
    }

    // compute challenge x
    std::string str1 = "";
    str1 += proof.A.ToByteString();
    str1 += proof.C.ToByteString();
    str1 += proof.D.ToByteString();
    BigInt x = Hash::StringToBigInt(str1);

    // prepare vec_x^m
    std::vector<BigInt> exp_x(m+1);
    exp_x[0] = bn_1;
    for(auto i=1;i<=m;i++){
        exp_x[i] = exp_x[i-1] * x % order;
    }

    // compute zG
    BigInt zG_leftleft = bn_0;
    for(auto i=0;i<k;i++){
        zG_leftleft = (zG_leftleft + vec_e_k[i] * witness.vec_r[i] % order) % order;
    }

    BigInt zG_right = bn_0;
    for(auto d=0;d<m;d++){
        zG_right = (zG_right + vec_rho[d] * exp_x[d] % order) % order;
    }

    proof.zG = (zG_leftleft * exp_x[m] % order - zG_right + order) % order;

    // compute zP   
    BigInt zP_right = bn_0;
    for(auto d=0;d<m;d++){
        zP_right = (zP_right + vec_alpha[d] * exp_x[d] % order) % order;
    }

    proof.zP = (rP * exp_x[m] % order + zP_right) % order;

    // Log Bit Proof
    // set pp
    LogBit::PP logbit_pp = LogBit::Setup(N);
    logbit_pp.g = pp.g;
    logbit_pp.h = pp.h;
    logbit_pp.u = pp.u;
    logbit_pp.u_new = pp.u_new;
    logbit_pp.vec_g = pp.vec_g;
    logbit_pp.vec_h = pp.vec_h;

    // set instance
    LogBit::Instance logbit_instance;
    logbit_instance.P = proof.P;

    // set witness
    LogBit::Witness logbit_witness;
    logbit_witness.r = rP;
    logbit_witness.vec_a = vec_s;
    logbit_witness.vec_b.resize(N);
    for(auto i=0;i<N;i++){
        logbit_witness.vec_b[i] = (vec_s[i] - bn_1 + order) % order;
    }

    // set transcript_str
    std::string logbit_transcript_str = "";

    // compute log_proof
    proof.logbit_proof = LogBit::Prove(logbit_pp,logbit_instance,logbit_witness,logbit_transcript_str);

    return proof;

}

bool Verify(PP &pp, Instance &instance, std::string &transcript_str, Proof &proof)
{
    size_t N = pp.vec_g.size();
    size_t m = log2(N);
    size_t k = pp.k;
    size_t mk = m * k;

    // computer the challenge e
    std::string str = "";
    str += proof.A.ToByteString();
    str += proof.B.ToByteString();
    str += proof.C.ToByteString();
    str += proof.D.ToByteString();
    str += proof.P.ToByteString();

    BigInt e = Hash::StringToBigInt(str);

    // compute challenge x
    std::string str1 = "";
    str1 += proof.A.ToByteString();
    str1 += proof.C.ToByteString();
    str1 += proof.D.ToByteString();
    BigInt x = Hash::StringToBigInt(str1);

    // compute p_i,j(x)
    std::vector<std::vector<BigInt>> P_ij;
    for(auto i=0;i<k;i++){
        //compute p_j(x)
        std::vector<BigInt> vec_P = GenRandomBigIntVectorLessThan(N,order); 

        for(auto j = 0; j < N; j++){
            vec_P[j] = bn_1;
            std::vector<size_t> vec_index = Decompose(j, 2, m); 
    
            for(auto b = 0; b < m; b++){        
                if(vec_index[b] == 1){
                    vec_P[j] = vec_P[j] * proof.vec_f[i*m+b] % order;
                }
                else{
                    vec_P[j] = vec_P[j] * ((x - proof.vec_f[i*m+b] + order) % order) % order;
                }    
            } 
        }
        P_ij.emplace_back(vec_P);
    }

    // prepare vec_x^m
    std::vector<BigInt> exp_x(m+1);
    exp_x[0] = bn_1;
    for(auto i=1;i<=m;i++){
        exp_x[i] = exp_x[i-1] * x % order;
    }
    
    // compute vec_e_k
    std::vector<BigInt> vec_e_k = GenBigIntPowerVector(k, e);

    // compute P_new
    ECPoint P_new_l = pp.vec_h[0];
    for(auto i=1;i<N;i++){
        P_new_l += pp.vec_h[i];
    } 
    P_new_l = P_new_l * (- exp_x[m] );

    ECPoint P_new_m = P_new_l;
    for(auto j=0;j<N;j++){
        BigInt index = bn_0;
        for(auto i=0;i<k;i++){
            index = (index + P_ij[i][j]) % order;
        }
        P_new_m += (pp.vec_g[j] + pp.vec_h[j]) * index;
    }

    ECPoint P_new_r = P_new_m;
    for(auto d=0;d<m;d++){
        P_new_r += proof.vec_C_p[d] * (- exp_x[d]);
    }

    ECPoint P_new = P_new_r;

    // compute G
    ECPoint G_l ;
    G_l.SetInfinity();
    for(auto j=0;j<N;j++){
        BigInt index = bn_0;
        for(auto i=0;i<k;i++){
            index = (index + vec_e_k[i] * P_ij[i][j] % order) % order;
        }
        G_l += instance.vec_c[j] * index;
    }

    ECPoint G_r = proof.vec_C_c[0] * (- exp_x[0]);
    for(auto d=1;d<m;d++){
        G_r += proof.vec_C_c[d] * (- exp_x[d]);
    }

    ECPoint G = G_l + G_r;

    // start verify
    std::vector<bool> vec_condition(4, false);

    // check condition 1
    LogBit::PP logbit_pp = LogBit::Setup(N);
    logbit_pp.g = pp.g;
    logbit_pp.h = pp.h;
    logbit_pp.u = pp.u;
    logbit_pp.u_new = pp.u_new;
    logbit_pp.vec_g = pp.vec_g;
    logbit_pp.vec_h = pp.vec_h;

    LogBit::Instance logbit_instance;
    logbit_instance.P = proof.P;

    std::string logbit_transcript_str = "";

    vec_condition[0] = LogBit::Verify(logbit_pp,logbit_instance,logbit_transcript_str,proof.logbit_proof);

    // check condition 2
    LinBit::PP linbit_pp = LinBit::Setup(mk);
    linbit_pp.h = pp.u;
    linbit_pp.vec_g = pp.vec_g_mk;
    
    LinBit::Instance linbit_instance; 
    linbit_instance.P = proof.B;
    
    LinBit::Proof linbit_proof;
    linbit_proof.A = proof.A;
    linbit_proof.C = proof.C;
    linbit_proof.D = proof.D;
    linbit_proof.zA = proof.zA;
    linbit_proof.zC = proof.zC;
    linbit_proof.vec_f.resize(mk);
    for(auto i=0;i<mk;i++){
        linbit_proof.vec_f[i] = proof.vec_f[i];
    }
    
    std::string transcript_str_linbit = "";
    vec_condition[1] = LinBit::Verify(linbit_pp, linbit_instance, transcript_str_linbit, linbit_proof);
    
    // check condition 3
    vec_condition[2] = (G == pp.g * proof.zG);

    // check condition 4
    vec_condition[3] = (proof.P * exp_x[m] + P_new.Invert() == pp.u * proof.zP);

    // check result
    bool Validity = vec_condition[0] && vec_condition[1] && vec_condition[2] && vec_condition[3] ;


    #ifdef DEBUG
    for(auto i = 0; i < 4; i++){
        std::cout << std::boolalpha << "Condition "<< std::to_string(i) <<" (Koon proof) = " 
                  << vec_condition[i] << std::endl; 
    }

    if (Validity){ 
        std::cout << "NIZK proof for Koon accepts >>>" << std::endl; 
    } else {
        std::cout << "NIZK proof for Koon rejects >>>" << std::endl; 
    }
    #endif

    return Validity;
}

}

#endif



