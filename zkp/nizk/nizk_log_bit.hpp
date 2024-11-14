
#ifndef NIZK_LOG_BIT_HPP_
#define NIZK_LOG_BIT_HPP_

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"
#include "../../pke/twisted_exponential_elgamal.hpp"
#include "../../commitment/pedersen.hpp"
#include "../../utility/polymul.hpp"
#include "../bulletproofs/innerproduct_proof.hpp"


namespace LogBit{

using Serialization::operator<<; 
using Serialization::operator>>; 

std::vector<size_t> Decompose(size_t l, size_t n, size_t m)
{
    std::vector<size_t> vec_index(m); 
    for(auto j = 0; j < m; j++){
        vec_index[j] = l % n;  
        l = l / n; 
    }
    return vec_index;  
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

struct PP
{
    std::vector<ECPoint> vec_g;
    std::vector<ECPoint> vec_h;
    ECPoint g,h,u;
    // for ip
    ECPoint u_new;
    
};

// structure of instance
struct Instance
{
    ECPoint P;
    // k in [0,n]
    size_t k;
};

// structure of witness 
struct Witness
{
    std::vector<BigInt> vec_a;
    std::vector<BigInt> vec_b;
    BigInt r;
};


// structure of proof 
struct Proof
{
    ECPoint A;
    ECPoint T1;
    ECPoint T2;
    BigInt t;
    BigInt tau;
    BigInt mu;
    InnerProduct::Proof ip_proof;
};
 
std::ofstream &operator<<(std::ofstream &fout, const Proof &proof)
{
    fout << proof.A << proof.T1 << proof.T2 << proof.t << proof.tau << proof.mu;
    return fout;  
}

std::ifstream &operator>>(std::ifstream &fin, Proof &proof)
{
    fin >> proof.A >> proof.T1 >> proof.T2 >> proof.t >> proof.tau >> proof.mu;
    return fin; 
}


/* Setup algorithm */ 
PP Setup(size_t N)
{
    PP pp; 
    
    pp.g = generator;
    pp.h = Hash::StringToECPoint(pp.g.ToByteString());
    pp.vec_g = GenRandomECPointVector(N);
    pp.vec_h = GenRandomECPointVector(N);
    pp.u = GenRandomGenerator();
    pp.u_new = GenRandomGenerator();

    return pp; 
}
  

Proof Prove(PP &pp, Instance &instance, Witness &witness, std::string &transcript_str)
{
    Proof proof;

    size_t N = pp.vec_g.size();

    // choose alpha beta rho
    std::vector<BigInt> vec_alpha = GenRandomBigIntVectorLessThan(N,order);
    std::vector<BigInt> vec_beta = GenRandomBigIntVectorLessThan(N,order);
    BigInt rho = GenRandomBigIntLessThan(order);

    // compute A
    proof.A = pp.u * rho;
    for(auto i=0;i<N;i++){
        proof.A += pp.vec_g[i] * vec_alpha[i];
        proof.A += pp.vec_h[i] * vec_beta[i];
    }

    // compute y,z
    transcript_str += proof.A.ToByteString();
    BigInt y = Hash::StringToBigInt(transcript_str);

    transcript_str += proof.A.ToByteString();
    BigInt z = Hash::StringToBigInt(transcript_str);

    // prepare the vector polynomials -- APGC Page 10 Figure 2
    
    // compute l(X)
    std::vector<BigInt> vec_z_unary(N, z);
    std::vector<BigInt> poly_ll0 = BigIntVectorModSub(witness.vec_a, vec_z_unary, BigInt(order));  
    std::vector<BigInt> poly_ll1(N); 
    poly_ll1 = vec_alpha; 

    // compute r(X)     
    std::vector<BigInt> vec_y_power = GenBigIntPowerVector(N, y); 
    std::vector<BigInt> vec_zz_temp = BigIntVectorModAdd(vec_z_unary, witness.vec_b, BigInt(order)); 
    std::vector<BigInt> poly_rr0 = BigIntVectorModProduct(vec_y_power, vec_zz_temp, BigInt(order)); 
    std::vector<BigInt> poly_rr1 = BigIntVectorModProduct(vec_y_power, vec_beta, BigInt(order));

    // compute t(X) 
    BigInt t0 = BigIntVectorModInnerProduct(poly_ll0, poly_rr0, BigInt(order)); 
    BigInt bn_temp1 = BigIntVectorModInnerProduct(poly_ll1, poly_rr0, BigInt(order)); 
    BigInt bn_temp2 = BigIntVectorModInnerProduct(poly_ll0, poly_rr1, BigInt(order));
    BigInt t1 = (bn_temp1 + bn_temp2) % BigInt(order);  
    BigInt t2 = BigIntVectorModInnerProduct(poly_ll1, poly_rr1, BigInt(order));

    // compute T1,T2
    BigInt tau1 = GenRandomBigIntLessThan(order);
    BigInt tau2 = GenRandomBigIntLessThan(order);

    proof.T1 = pp.g * t1 + pp.h * tau1; 
    proof.T2 = pp.g * t2 + pp.h * tau2; 

    // compute x
    transcript_str += proof.T1.ToByteString();
    transcript_str += proof.T2.ToByteString();           

    BigInt x = Hash::StringToBigInt(transcript_str);

    BigInt x_square = x.ModSquare(order);   

    // compute the value of l(x) and r(x) at point x
    vec_zz_temp = BigIntVectorModScalar(poly_ll1, x, BigInt(order));
    std::vector<BigInt> llx = BigIntVectorModAdd(poly_ll0, vec_zz_temp, BigInt(order));

    vec_zz_temp = BigIntVectorModScalar(poly_rr1, x, BigInt(order)); 
    std::vector<BigInt> rrx = BigIntVectorModAdd(poly_rr0, vec_zz_temp, BigInt(order));

    // compute t
    proof.t = BigIntVectorModInnerProduct(llx, rrx, BigInt(order)); 
 
    // compute taux
    proof.tau = ((tau1 * x % order) + (tau2 * x_square % order)) % order;

    // compute mu 
    proof.mu = (witness.r + rho * x % order) % order;

    // compute vec_h_new
    BigInt y_inverse = y.ModInverse(order);
    std::vector<BigInt> vec_y_inverse_power = GenBigIntPowerVector(N, y_inverse);

    std::vector<ECPoint> vec_h_new(N);
    for(auto i=0;i<N;i++){
        vec_h_new[i] = pp.vec_h[i] * vec_y_inverse_power[i];
    }
    
    // compute B
    ECPoint B = instance.P + proof.A * x;
    for(auto i=0;i<N;i++){
        B += vec_h_new[i] * (z * vec_y_power[i] % order);
        B += pp.vec_g[i] * z.ModNegate(order);
    }
    
    // prepare u for ip
    ECPoint P_temp = B + (pp.u * proof.mu).Invert();

    transcript_str += P_temp.ToByteString() + proof.t.ToByteString();
    BigInt x_new = Hash::StringToBigInt(transcript_str);

    ECPoint u_ip = pp.u_new * x_new;

    // inner product proof
    // set ip pp
    InnerProduct::PP ip_pp = InnerProduct::Setup(N, true); 
    ip_pp.vec_g = pp.vec_g;
    ip_pp.vec_h = vec_h_new;
    ip_pp.u = u_ip;

    // set ip witness
    InnerProduct::Witness ip_witness;
    ip_witness.vec_a = llx; 
    ip_witness.vec_b = rrx;

    // set ip instance    
    InnerProduct::Instance ip_instance;
    ip_instance.P = P_temp + ip_pp.u * proof.t;  

    InnerProduct::Prove(ip_pp, ip_instance, ip_witness, transcript_str, proof.ip_proof);  
    
    return proof;

}


bool Verify(PP &pp, Instance &instance, std::string &transcript_str, Proof &proof)
{
    size_t N = pp.vec_g.size();
    
    std::vector<bool> vec_condition(2, false);

    // compute y,z,x
    transcript_str += proof.A.ToByteString();
    BigInt y = Hash::StringToBigInt(transcript_str);

    transcript_str += proof.A.ToByteString();
    BigInt z = Hash::StringToBigInt(transcript_str);
    
    transcript_str += proof.T1.ToByteString() + proof.T2.ToByteString();
    BigInt x = Hash::StringToBigInt(transcript_str);

    // compute y power
    std::vector<BigInt> vec_y_power = GenBigIntPowerVector(N, y); // y^nm

    // compute vec_h_new
    BigInt y_inverse = y.ModInverse(order);
    std::vector<BigInt> vec_y_inverse_power = GenBigIntPowerVector(N, y_inverse);

    std::vector<ECPoint> vec_h_new(N);
    for(auto i=0;i<N;i++){
        vec_h_new[i] = pp.vec_h[i] * vec_y_inverse_power[i];
    }   

    // compute B
    ECPoint B = instance.P + proof.A * x;
    for(auto i=0;i<N;i++){
        B += vec_h_new[i] * (z * vec_y_power[i] % order);
        B += pp.vec_g[i] * z.ModNegate(order);
    }

    // prepare u for ip
    ECPoint P_temp = B + (pp.u * proof.mu).Invert();

    transcript_str += P_temp.ToByteString() + proof.t.ToByteString();
    BigInt x_new = Hash::StringToBigInt(transcript_str);

    ECPoint u_ip = pp.u_new * x_new;
    
    // prepare
    std::vector<BigInt> vec_1_power(N, bn_1); 
    BigInt z_square = z.ModSquare(order);

    // compute delta_yz -- APGC Page 8 eq(2): 2^n -> 0 for vec_w = 0    
    BigInt bn_temp1 = BigIntVectorModInnerProduct(vec_1_power, vec_y_power, BigInt(order)); 
    BigInt delta_yz = z.ModSub(z_square, order) * bn_temp1; 

    // check condition 1
    
    // inner product proof
    // set ip pp
    InnerProduct::PP ip_pp = InnerProduct::Setup(N, true);
    ip_pp.vec_g = pp.vec_g;
    ip_pp.vec_h = vec_h_new;
    ip_pp.u = u_ip;

    // set ip instance
    InnerProduct::Instance ip_instance;
    ip_instance.P = P_temp + ip_pp.u * proof.t;    

 
    vec_condition[0] = InnerProduct::Verify(ip_pp, ip_instance, transcript_str, proof.ip_proof);  
    
    // check condition 2
    // left
    ECPoint left = pp.g * proof.t + pp.h * proof.tau;

    // right
    ECPoint right = pp.g * delta_yz + proof.T1 * x + proof.T2 * x.ModSquare(order);

    vec_condition[1] = (left == right);

    bool Validity = vec_condition[0] && vec_condition[1];

    #ifdef DEBUG
    for(auto i = 0; i < 2; i++){
        std::cout << std::boolalpha << "Condition "<< std::to_string(i) <<" (Log Bit proof) = " 
                  << vec_condition[i] << std::endl; 
    }

    if (Validity){ 
        std::cout << "NIZK proof for Log Bit accepts >>>" << std::endl; 
    } else {
        std::cout << "NIZK proof for Log Bit rejects >>>" << std::endl; 
    }
    #endif

    return Validity;
}

}

#endif



