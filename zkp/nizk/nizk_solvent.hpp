#ifndef NIZK_SOLVENT_HPP_
#define NIZK_SOLVENT_HPP_

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"
#include "../../pke/twisted_exponential_elgamal.hpp"
#include "../../commitment/pedersen.hpp"
#include "../../utility/polymul.hpp"
#include "../../zkp/bulletproofs/bullet_proof.hpp"
#include "../../zkp/nizk/nizk_ACF_1oon_range_for_solvent.hpp"

namespace Solvent{

using Serialization::operator<<; 
using Serialization::operator>>; 

struct PP
{
    size_t N;
    std::vector<ECPoint> vec_g_linbit;
    std::vector<ECPoint> vec_g_range1;
    std::vector<ECPoint> vec_g_range2;
    std::vector<ECPoint> vec_h1;
    std::vector<ECPoint> vec_h2;
    ECPoint g,h,u,u_new;
    std::vector<ECPoint> vec_g_acf;
    std::vector<ECPoint> vec_h_acf;
    std::vector<ECPoint> vec_u_acf;
    std::vector<ECPoint> vec_ip_g;
    std::vector<ECPoint> vec_ip_h;
    
};

// structure of instance
struct Instance
{
    std::vector<ECPoint> vec_C;
    std::vector<ECPoint> vec_pk;
    std::vector<ECPoint> vec_CL;
    std::vector<ECPoint> vec_CR;
    std::vector<ECPoint> vec_C_L;
    std::vector<ECPoint> vec_C_R;
};

// structure of witness 
struct Witness
{
    BigInt sk;
    size_t l;
    BigInt v;
    BigInt vl;
    BigInt rl;
};

// structure of proof 
struct Proof
{
    ECPoint A,B,C,D;
    std::vector<ECPoint> vec_Cc;
    std::vector<ECPoint> vec_Cs;
    ECPoint V1,V2,S1,S2;
    ECPoint T11,T12,T21,T22;
    std::vector<BigInt> vec_f;
    BigInt zA,zC,zG,zS;
    BigInt mu1,mu2,t1,t2,tau1,tau2;
    std::vector<ECPoint> vec_L1, vec_L2, vec_R1, vec_R2;
    BigInt a1,a2,b1,b2;
    ACF_1oon_range::Proof ACF_proof; 
};
 
std::ofstream &operator<<(std::ofstream &fout, const Proof &proof)
{
    fout << proof.A << proof.B << proof.C << proof.D ;
    fout << proof.vec_Cc << proof.vec_Cs;
    fout << proof.V1 << proof.V2 << proof.S1 << proof.S2 ;
    fout << proof.T11 << proof.T12 << proof.T21 << proof.T22 ;
    fout << proof.vec_f;
    fout << proof.zA << proof.zC << proof.zG << proof.zS;
    fout << proof.mu1 << proof.mu2 << proof.t1 << proof.t2 << proof.tau1 << proof.tau2 ;
    fout << proof.vec_L1 << proof.vec_L2 << proof.vec_R1 << proof.vec_R2 ;
    fout << proof.a1 << proof.a2 << proof.b1 << proof.b2 ;
    return fout;
}

std::ifstream &operator>>(std::ifstream &fin, Proof &proof)
{
    fin >> proof.A >> proof.B >> proof.C >> proof.D ;
    fin >> proof.vec_Cc >> proof.vec_Cs;
    fin >> proof.V1 >> proof.V2 >> proof.S1 >> proof.S2 ;
    fin >> proof.T11 >> proof.T12 >> proof.T21 >> proof.T22 ;
    fin >> proof.vec_f;
    fin >> proof.zA >> proof.zC >> proof.zG >> proof.zS ;
    fin >> proof.mu1 >> proof.mu2 >> proof.t1 >> proof.t2 >> proof.tau1 >> proof.tau2 ;
    fin >> proof.vec_L1 >> proof.vec_L2 >> proof.vec_R1 >> proof.vec_R2 ;
    fin >> proof.a1 >> proof.a2 >> proof.b1 >> proof.b2 ;
    return fin;
}

/* Setup algorithm */ 
PP Setup(size_t N)
{
    PP pp; 
    pp.N = N;
    pp.vec_g_linbit = GenRandomECPointVector(log2(N));
    pp.vec_g_range1 = GenRandomECPointVector(32);
    pp.vec_g_range2 = GenRandomECPointVector(32);
    pp.vec_h1 = GenRandomECPointVector(32);
    pp.vec_h2 = GenRandomECPointVector(32);
    pp.g = GenRandomGenerator(); 
    pp.h = GenRandomGenerator(); 
    pp.u = GenRandomGenerator(); 
    pp.u_new = GenRandomGenerator(); 
    pp.vec_g_acf = GenRandomECPointVector(N);
    pp.vec_h_acf = GenRandomECPointVector(N);
    pp.vec_u_acf = GenRandomECPointVector(3*N);
    pp.vec_ip_h = GenRandomECPointVector(32);
    pp.vec_ip_g = GenRandomECPointVector(32);

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

Proof Prove(PP &pp, Instance &instance, Witness &witness, std::string &transcript_str)
{
    Proof proof;

    size_t N = pp.N;
    size_t m = log2(N);

    ECPoint init;
    init.SetInfinity();

    // run ACF 1oon dgh range
    ACF_1oon_range::PP pp_acf = ACF_1oon_range::Setup(N);
    pp_acf.g = pp.g;
    pp_acf.h = pp.h;
    pp_acf.u = pp.u;
    pp_acf.u_new = pp.u_new;
    pp_acf.vec_g = pp.vec_g_acf;
    pp_acf.vec_h = pp.vec_h_acf;
    pp_acf.vec_u = pp.vec_u_acf;
    pp_acf.vec_ip_g = pp.vec_ip_g;
    pp_acf.vec_ip_h = pp.vec_ip_h;

    ACF_1oon_range::Instance instance_acf;
    instance_acf.vec_pk = instance.vec_pk;
    instance_acf.vec_CL = instance.vec_CL;
    instance_acf.vec_CR = instance.vec_CR;
    instance_acf.vec_C_L = instance.vec_C_L;
    instance_acf.vec_C_R = instance.vec_C_R;

    ACF_1oon_range::Witness witness_acf;
    witness_acf.l = witness.l;
    witness_acf.vec_x.resize(2);
    witness_acf.vec_x[0] = witness.sk.ModInverse(order);
    witness_acf.vec_x[1] = witness.v;

    proof.ACF_proof = ACF_1oon_range::Prove(pp_acf,instance_acf,witness_acf,transcript_str);

    // decompose l
    std::vector<size_t> l_index = Decompose(witness.l,2,m);

    std::vector<BigInt> vec_l(m,bn_0);
    for(auto i=0;i<m;i++){
        if(l_index[i] == 1){
            vec_l[i] = bn_1;
        }
    } 

    // pick rB
    BigInt rB = GenRandomBigIntLessThan(order);

    // compute B
    proof.B = ECPointVectorMul(pp.vec_g_linbit, vec_l) + pp.u * rB;

    // run lin bit
    // pick rA,rC,rD,vec_alpha
    BigInt rA = GenRandomBigIntLessThan(order);
    BigInt rC = GenRandomBigIntLessThan(order);
    BigInt rD = GenRandomBigIntLessThan(order);
    std::vector<BigInt> vec_alpha = GenRandomBigIntVectorLessThan(m,order);

    // prepare for compute C,D
    std::vector<BigInt> vec_C_index(m,bn_0);
    std::vector<BigInt> vec_D_index(m,bn_0);
    for(auto i=0;i<m;i++){
        vec_C_index[i] = vec_alpha[i] * (bn_1 - bn_2 * l_index[i]) % order;
        vec_D_index[i] = (bn_0 - vec_alpha[i]) * vec_alpha[i] % order;
    }

    // compute A,C,D
    proof.A = ECPointVectorMul(pp.vec_g_linbit, vec_alpha) + pp.u * rA;
    proof.C = ECPointVectorMul(pp.vec_g_linbit, vec_C_index) + pp.u * rC;
    proof.D = ECPointVectorMul(pp.vec_g_linbit, vec_D_index) + pp.u * rD;

    // compute challenge
    transcript_str += proof.A.ToByteString() + proof.C.ToByteString() + proof.D.ToByteString();
    BigInt x0 = Hash::StringToBigInt(transcript_str);

    // compute vec_f,zA,zC
    std::vector<BigInt> vec_temp = BigIntVectorModScalar(vec_l, x0,order);
    proof.vec_f = BigIntVectorModAdd(vec_temp, vec_alpha, order);
    proof.zA = (rB * x0 + rA) % order;
    proof.zC = (rC * x0 + rD) % order;

    //compute p_j_d
    std::vector<std::vector<BigInt>> P; 
    for(auto j = 0; j < N; j++){
        std::vector<std::vector<BigInt>> A(m, std::vector<BigInt>(2));        
        // prepare m ploynomial of form ax+b
        std::vector<size_t> vec_index = Decompose(j, 2, m); 
        for(auto b = 0; b < m; b++){      
            if(vec_index[b] == 1){
                A[b][0] = vec_alpha[b];
                A[b][1] = l_index[b];
            }
            else{
                A[b][0] = bn_0 - vec_alpha[b];
                A[b][1] = bn_1 - l_index[b];
            }    
        } 
        std::vector<BigInt> p_j = PolyMul(A);    
        P.emplace_back(p_j); 
    }

    // pick vec_alpha_d,vec_beta_d
    std::vector<BigInt> vec_alpha_d = GenRandomBigIntVectorLessThan(m,order);
    std::vector<BigInt> vec_beta_d = GenRandomBigIntVectorLessThan(m,order);

    // compute vec_Cc,vec_Cs
    proof.vec_Cc.resize(m);
    proof.vec_Cs.resize(m);
    for(auto d=0;d<m;d++){
        proof.vec_Cc[d] = pp.g * vec_alpha_d[d];
        proof.vec_Cs[d] = pp.u * vec_beta_d[d];
        for(auto j=0;j<N;j++){
            proof.vec_Cc[d] += instance.vec_CR[j] * P[j][d];
            proof.vec_Cs[d] += (pp.vec_g_acf[j] + pp.vec_h_acf[j]) * P[j][d];
        }
    }

    // compute vec_x0_power
    std::vector<BigInt> vec_x0_power = GenBigIntPowerVector(m+1,x0);

    // compute zG，zS
    proof.zG = witness.rl * vec_x0_power[m] % order;
    proof.zS = proof.ACF_proof.rS * vec_x0_power[m] % order;
    for(auto d=0;d<m;d++){
        proof.zG = (proof.zG - vec_alpha_d[d] * vec_x0_power[d] + order) % order;
        proof.zS = (proof.zS + vec_beta_d[d] * vec_x0_power[d]) % order;
    }
    
    // compute tau
    BigInt tau = (proof.zG * vec_x0_power[m].ModInverse(order) % order).ModNegate(order);
    
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
    First_Bullet_witness.v[0] = witness.vl.ModNegate(order);

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
    Second_Bullet_witness.v[0] = witness.vl.ModNegate(order) - bn_1;

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

    #ifdef DEBUG
    std::cout << "Sdr Trans Proof Generation Finishes >>>" << std::endl;
    #endif 

    return proof;

}

bool Verify(PP &pp, Instance &instance, std::string &transcript_str, Proof &proof)
{
    size_t N = pp.N;
    size_t m = log2(N);
    ECPoint init;
    init.SetInfinity();

    std::vector<bool> vec_condition(5, false);

    // check ACF 1oon dgh range
    ACF_1oon_range::PP pp_acf = ACF_1oon_range::Setup(N);
    pp_acf.g = pp.g;
    pp_acf.h = pp.h;
    pp_acf.u = pp.u;
    pp_acf.u_new = pp.u_new;
    pp_acf.vec_g = pp.vec_g_acf;
    pp_acf.vec_h = pp.vec_h_acf;
    pp_acf.vec_u = pp.vec_u_acf;
    pp_acf.vec_ip_g = pp.vec_ip_g;
    pp_acf.vec_ip_h = pp.vec_ip_h;

    ACF_1oon_range::Instance instance_acf;
    instance_acf.vec_pk = instance.vec_pk;
    instance_acf.vec_CL = instance.vec_CL;
    instance_acf.vec_CR = instance.vec_CR;
    instance_acf.vec_C_L = instance.vec_C_L;
    instance_acf.vec_C_R = instance.vec_C_R;

    vec_condition[0] = ACF_1oon_range::Verify(pp_acf,instance_acf,transcript_str,proof.ACF_proof);

    // recover x0
    transcript_str += proof.A.ToByteString() + proof.C.ToByteString() + proof.D.ToByteString();
    BigInt x0 = Hash::StringToBigInt(transcript_str);

    // compute vec_x0_power
    std::vector<BigInt> vec_x0_power = GenBigIntPowerVector(m+1,x0);

    // compute x0m_inv
    BigInt x0m_inv = vec_x0_power[m].ModInverse(order);

    //compute p_j(x)
    std::vector<BigInt> vec_P(N,bn_1); 
    for(auto j = 0; j < N; j++){
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

    // prepare for S1
    ECPoint S11 = init;
    for(auto j=0;j<N;j++){
        S11 += pp.vec_h_acf[j];
    }
    S11 = S11 * (-vec_x0_power[m]);

    // compute G1,S1
    ECPoint G1 = init;
    ECPoint S1 = S11;
    for(auto j=0;j<N;j++){
        G1 += instance.vec_CR[j] * vec_P[j];
        S1 += (pp.vec_g_acf[j] + pp.vec_h_acf[j]) * vec_P[j];
    }
    for(auto d=0;d<m;d++){
        G1 += proof.vec_Cc[d] * (bn_0 - vec_x0_power[d] + order);
        S1 += proof.vec_Cs[d] * (bn_0 - vec_x0_power[d] + order);
    }

    // check lin bit
    std::vector<BigInt> vec_of_f(m,bn_0);
    for(auto i=0;i<m;i++){
        vec_of_f[i] = proof.vec_f[i] * (x0 - proof.vec_f[i]) % order;
    }
    ECPoint R1,R2;
    R1 = ECPointVectorMul(pp.vec_g_linbit, proof.vec_f);
    R2 = ECPointVectorMul(pp.vec_g_linbit, vec_of_f);

    bool LB1 = ((proof.B * x0 + proof.A) == (R1 + pp.u * proof.zA));
    bool LB2 = ((proof.C * x0 + proof.D) == (R2 + pp.u * proof.zC));

    vec_condition[1] = LB1 && LB2;

    // check S,S1
    ECPoint LEFT = proof.ACF_proof.S * vec_x0_power[m] + S1.Invert();
    ECPoint RIGHT = pp.u * proof.zS;

    vec_condition[2] = (LEFT == RIGHT);

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
    First_Bullet_instance.C[0] = G1 * x0m_inv.ModNegate(order);
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

    vec_condition[3] = Bullet::Verify(First_Bullet_pp,First_Bullet_instance,First_Bullet_str,First_Bullet_proof);

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
    Second_Bullet_instance.C[0] = G1 * x0m_inv.ModNegate(order) + pp.h.Invert();

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

    vec_condition[4] = Bullet::Verify(Second_Bullet_pp,Second_Bullet_instance,Second_Bullet_str,Second_Bullet_proof);

    bool Validity = vec_condition[0] && vec_condition[1] && vec_condition[2] && vec_condition[3] && vec_condition[4];

    #ifdef DEBUG
    for(auto i = 0; i < 5; i++){
        std::cout << std::boolalpha << "Condition "<< std::to_string(i) <<" (Solvent proof) = " 
                  << vec_condition[i] << std::endl; 
    }

    if (Validity){ 
        std::cout << "NIZK proof for Solvent accepts >>>" << std::endl; 
    } else {
        std::cout << "NIZK proof for Solvent rejects >>>" << std::endl; 
    }
    #endif

    return Validity;
}

}

#endif





