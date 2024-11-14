
#ifndef BULLET_PROOF2_HPP_
#define BULLET_PROOF2_HPP_

#include "innerproduct_proof.hpp" 

namespace Bullet{

using Serialization::operator<<; 
using Serialization::operator>>; 

// define the structure of Bulletproofs
struct PP
{
    size_t RANGE_LEN; 
    size_t LOG_RANGE_LEN; 
    size_t MAX_AGG_NUM; // number of sub-argument (for now, we require m to be the power of 2)

    ECPoint g, h;
    ECPoint u; // used for inside innerproduct statement

    BigInt e;
    std::vector<ECPoint> vec_g, vec_h; // the pp of innerproduct part    
};

struct Instance
{
    std::vector<ECPoint> C;  // Ci = g^ri h^vi: length = AGG_NUM
}; 

struct Witness
{
    std::vector<BigInt> r; // length = AGG_NUM
    std::vector<BigInt> v; 
}; 

struct Proof
{
    ECPoint A, S, T1, T2;  
    BigInt taux, mu, tx; 
    InnerProduct::Proof ip_proof;    
};

std::ofstream &operator<<(std::ofstream &fout, const Proof &proof)
{
    fout << proof.A << proof.S << proof.T1 << proof.T2 << proof.taux << proof.mu << proof.tx; 
    fout << proof.ip_proof; 
    return fout; 
}

std::ifstream &operator>>(std::ifstream &fin, InnerProduct::Proof &proof)
{
    fin >> proof.A >> proof.S >> proof.T1 >> proof.T2 >> proof.taux >> proof.mu >> proof.tx; 
    fin >> proof.ip_proof; 
    return fin; 
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

PP Setup(size_t &RANGE_LEN, size_t &MAX_AGG_NUM)
{
    PP pp; 
    pp.RANGE_LEN = RANGE_LEN; 
    pp.LOG_RANGE_LEN = size_t(log2(RANGE_LEN)); 
    pp.MAX_AGG_NUM = MAX_AGG_NUM; 
 
    pp.g = generator; 
    pp.h = Hash::StringToECPoint(pp.g.ToByteString()); 
    pp.u = GenRandomGenerator();
    pp.e = GenRandomBigIntLessThan(order);
    pp.vec_g = GenRandomECPointVector(RANGE_LEN*MAX_AGG_NUM);
    pp.vec_h = GenRandomECPointVector(RANGE_LEN*MAX_AGG_NUM);

    return pp; 
}

// statement C = g^r h^v and v \in [0, 2^n-1]
void Prove(PP &pp, Instance &instance, Witness &witness, std::string &transcript_str, Proof &proof)
{ 
    size_t n = pp.MAX_AGG_NUM;
    size_t LEN = pp.RANGE_LEN * n; // LEN = mn

    std::vector<BigInt> vec_aL(LEN);  
    std::vector<BigInt> vec_aR(LEN);
 
    std::vector<BigInt> vec_1_power(LEN, bn_1); // vec_unary = 1^nm

    for (auto i = 0; i < n; i++)
    {
        for(auto j = 0; j < pp.RANGE_LEN; j++)
        {
            if(witness.v[i].GetTheNthBit(j) == 1){
                vec_aL[i * pp.RANGE_LEN + j] = bn_1;  
            }
            else{
                vec_aL[i * pp.RANGE_LEN + j] = bn_0; 
            } 
        }
    }

    vec_aR = BigIntVectorModSub(vec_aL, vec_1_power,  BigInt(order)); // Eq (42) -- aR = aL - 1^n

    // prepare vec_A and vec_a for multi-exponention (used hereafter)
    
    // Eq (44) -- compute A = H^alpha g^aL h^aR (commitment to aL and aR)
    BigInt alpha = GenRandomBigIntLessThan(order); 

    std::vector<ECPoint> vec_A(2*LEN+1); 
    std::copy(pp.vec_g.begin(), pp.vec_g.begin()+LEN, vec_A.begin()); 
    std::copy(pp.vec_h.begin(), pp.vec_h.begin()+LEN, vec_A.begin()+LEN); 
    vec_A[2*LEN] = pp.h; 

    std::vector<BigInt> vec_a(2*LEN+1); 
    std::copy(vec_aL.begin(), vec_aL.begin()+LEN, vec_a.begin()); 
    std::copy(vec_aR.begin(), vec_aR.begin()+LEN, vec_a.begin()+LEN); 
    vec_a[2*LEN] = alpha;

    proof.A = ECPointVectorMul(vec_A, vec_a); // Eq (44) 

    // pick sL, sR from Z_p^n (choose blinding vectors sL, sR)
    std::vector<BigInt> vec_sL = GenRandomBigIntVectorLessThan(LEN, order); 
    std::vector<BigInt> vec_sR = GenRandomBigIntVectorLessThan(LEN, order); 
    
    // Eq (47) compute S = H^alpha g^aL h^aR (commitment to sL and sR)
    BigInt rho = GenRandomBigIntLessThan(order); 

    std::copy(vec_sL.begin(), vec_sL.end(), vec_a.begin()); 
    std::copy(vec_sR.begin(), vec_sR.end(), vec_a.begin()+LEN); 
    vec_a[2*LEN] = rho; 

    proof.S = ECPointVectorMul(vec_A, vec_a); // Eq (47) 

    // Eq (49, 50) compute y and z
    transcript_str += proof.A.ToByteString(); 
    BigInt y = Hash::StringToBigInt(transcript_str);

    BigInt y_inverse = y.ModInverse(order);
     
    std::vector<BigInt> vec_y_inverse_power = GenBigIntPowerVector(LEN, y_inverse); // y^{-i+1}

    transcript_str += proof.S.ToByteString(); 
    BigInt z = Hash::StringToBigInt(transcript_str);

    BigInt z_square = z.ModSquare(order);
    BigInt z_cubic = (z * z_square) % order;
    
    std::vector<BigInt> vec_adjust_z_power(n+1); // generate z^{j+1} j \in [n] 
    vec_adjust_z_power[0] = z; 
    for (auto j = 1; j <= n; j++)
    {
        vec_adjust_z_power[j] = (z * vec_adjust_z_power[j-1]) % order; //pow(z, j+1, q); description below Eq (71)
    }  

    // prepare the vector polynomials
    
    // compute l(X) Eq (70)
    std::vector<BigInt> vec_z_unary(LEN, z); // z \cdot 1^nm

    std::vector<BigInt> poly_ll0 = BigIntVectorModSub(vec_aL, vec_z_unary, BigInt(order));  
    std::vector<BigInt> poly_ll1(LEN); 
    poly_ll1.assign(vec_sL.begin(), vec_sL.end()); 

    // compute r(X)     
    std::vector<BigInt> vec_y_power = GenBigIntPowerVector(LEN, y); // y^nm
    std::vector<BigInt> vec_zz_temp = BigIntVectorModAdd(vec_z_unary, vec_aR, BigInt(order)); // vec_t = aR + z1^nm
    std::vector<BigInt> poly_rr0 = BigIntVectorModProduct(vec_y_power, vec_zz_temp, BigInt(order)); // y^nm(aR + z1^nm)
    
    std::vector<BigInt> vec_short_2_power = GenBigIntPowerVector(pp.RANGE_LEN, bn_2); // 2^n
    std::vector<BigInt> vec_e_power = GenBigIntPowerVector(n,pp.e);

    for(auto i=0;i<n;i++){
        for(auto j=0;j<pp.RANGE_LEN;j++){
            vec_zz_temp[i*pp.RANGE_LEN+j] = (vec_e_power[i] * vec_short_2_power[j]) * z_square % order;
        }
    }
    
    poly_rr0 = BigIntVectorModAdd(poly_rr0, vec_zz_temp,BigInt(order));

    std::vector<BigInt> poly_rr1 = BigIntVectorModProduct(vec_y_power, vec_sR, BigInt(order)); //y^nsR X

    // compute t(X) 
    BigInt t0 = BigIntVectorModInnerProduct(poly_ll0, poly_rr0, BigInt(order)); 
    BigInt bn_temp1 = BigIntVectorModInnerProduct(poly_ll1, poly_rr0, BigInt(order)); 
    BigInt bn_temp2 = BigIntVectorModInnerProduct(poly_ll0, poly_rr1, BigInt(order));
    BigInt t1 = (bn_temp1 + bn_temp2) % BigInt(order);  
  
    BigInt t2 = BigIntVectorModInnerProduct(poly_ll1, poly_rr1, BigInt(order)); 

    // Eq (53) -- commit to t1, t2
    // P picks tau1 and tau2
    BigInt tau1 = GenRandomBigIntLessThan(order); 
    BigInt tau2 = GenRandomBigIntLessThan(order); 

    vec_A.clear(); vec_A = {pp.g, pp.h};
    // vec_A.clear(); vec_A = {pp.h, pp.g};
    
    vec_a.clear(); vec_a = {tau1, t1};  
    proof.T1 = ECPointVectorMul(vec_A, vec_a); //pp.g * tau1 + pp.h * t1; mul(tau1, pp.g, t1, pp.h);
    
    vec_a.clear(); vec_a = {tau2, t2};  
    proof.T2 = ECPointVectorMul(vec_A, vec_a); //pp.g * tau2 + pp.h * t2; mul(tau2, pp.g, t2, pp.h);    

    // Eq (56) -- compute the challenge x
    transcript_str += proof.T1.ToByteString() + proof.T2.ToByteString(); 
    BigInt x = Hash::StringToBigInt(transcript_str); 

    BigInt x_square = x.ModSquare(order);   

    // compute the value of l(x) and r(x) at point x
    vec_zz_temp = BigIntVectorModScalar(poly_ll1, x, BigInt(order));
    std::vector<BigInt> llx = BigIntVectorModAdd(poly_ll0, vec_zz_temp, BigInt(order));

    vec_zz_temp = BigIntVectorModScalar(poly_rr1, x, BigInt(order)); 
    std::vector<BigInt> rrx = BigIntVectorModAdd(poly_rr0, vec_zz_temp, BigInt(order)); 

    proof.tx = BigIntVectorModInnerProduct(llx, rrx, BigInt(order));  // Eq (60)  

    // compute taux
    proof.taux = (tau1 * x + tau2 * x_square + z_square * witness.r[0]) % order; //proof.taux = tau2*x_square + tau1*x; 

    // compute proof.mu = (alpha + rho*x) %q;  Eq (62)
    proof.mu = (alpha + rho * x) % order; 
    
    // transmit llx and rrx via inner product proof

    InnerProduct::PP ip_pp = InnerProduct::Setup(LEN, false); 
    ip_pp.vec_g.resize(LEN); 
    std::copy(pp.vec_g.begin(), pp.vec_g.begin()+LEN, ip_pp.vec_g.begin()); // ip_pp.vec_g = pp.vec_g

    ip_pp.vec_h.resize(LEN); 
    std::copy(pp.vec_h.begin(), pp.vec_h.begin()+LEN, ip_pp.vec_h.begin()); 
    ip_pp.vec_h = ECPointVectorProduct(ip_pp.vec_h, vec_y_inverse_power);  // ip_pp.vec_h = vec_h_new  

    transcript_str += x.ToByteString();  
    BigInt e = Hash::StringToBigInt(transcript_str);   

    InnerProduct::Witness ip_witness;
    ip_witness.vec_a = llx; // ip_witness.vec_a = llx
    ip_witness.vec_b = rrx; // ip_witness.vec_b = rrx

    InnerProduct::Instance ip_instance;
    ip_pp.u = pp.u * e; //ip_pp.u = u^e 

    vec_A.resize(2*LEN+1); 
    std::copy(ip_pp.vec_g.begin(), ip_pp.vec_g.end(), vec_A.begin()); 
    std::copy(ip_pp.vec_h.begin(), ip_pp.vec_h.end(), vec_A.begin()+LEN); 
    vec_A[2*LEN] = ip_pp.u;

    vec_a.resize(2*LEN+1); 
    std::copy(ip_witness.vec_a.begin(), ip_witness.vec_a.end(), vec_a.begin()); 
    std::copy(ip_witness.vec_b.begin(), ip_witness.vec_b.end(), vec_a.begin()+LEN); 
    vec_a[2*LEN] = proof.tx; 

    ip_instance.P = ECPointVectorMul(vec_A, vec_a);  
 
    InnerProduct::Prove(ip_pp, ip_instance, ip_witness, transcript_str, proof.ip_proof); 

    #ifdef DEBUG
        std::cout << "Bullet Proof Generation Succeeds >>>" << std::endl; 
    #endif
}

bool Verify(PP &pp, Instance &instance, std::string &transcript_str, Proof &proof)
{
    #ifdef DEBUG
        std::cout << "begin to check the proof" << std::endl; 
    #endif

    bool V1, V2, Validity; // variables for checking results

    transcript_str += proof.A.ToByteString(); 
    BigInt y = Hash::StringToBigInt(transcript_str);  //recover the challenge y
    BigInt y_inverse = y.ModInverse(order);  
    
    transcript_str += proof.S.ToByteString(); 
    BigInt z = Hash::StringToBigInt(transcript_str); // recover the challenge z

    BigInt z_minus = z.ModNegate(order); 
    BigInt z_square = z.ModSquare(order); // (z*z)%q; 
    BigInt z_cubic = (z * z_square) % order; 

    transcript_str += proof.T1.ToByteString() + proof.T2.ToByteString(); 
    BigInt x = Hash::StringToBigInt(transcript_str); 
    BigInt x_square = x.ModSquare(order);  // (x*x)%q;  //recover the challenge x from PI

    transcript_str += x.ToByteString(); 
    BigInt e = Hash::StringToBigInt(transcript_str);  // play the role of x_u

    // size_t n = instance.C.size();
    size_t n = pp.MAX_AGG_NUM;
    size_t LEN = pp.RANGE_LEN * n; // l = nm 
    std::vector<BigInt> vec_1_power(LEN, bn_1); // vec_unary = 1^nm
    std::vector<BigInt> vec_short_1_power(pp.RANGE_LEN, bn_1); 
    std::vector<BigInt> vec_2_power = GenBigIntPowerVector(LEN, bn_2);
    std::vector<BigInt> vec_short_2_power = GenBigIntPowerVector(pp.RANGE_LEN, bn_2);  
    std::vector<BigInt> vec_y_power = GenBigIntPowerVector(LEN, y); 

    std::vector<BigInt> vec_adjust_z_power(n+1); // generate z^{j+2} j \in [n]
    vec_adjust_z_power[0] = z; 
    for (auto j = 1; j <= n; j++)
    {
        vec_adjust_z_power[j] = (z * vec_adjust_z_power[j-1]) % order; 
    }  

    // compute sum_{j=1^m} z^{j+2}
    BigInt sum_z = bn_0; 
    for (auto j = 1; j <= n; j++)
    {
        sum_z += vec_adjust_z_power[j]; 
    }
    sum_z = (sum_z * z) % order;  

    // compute delta_yz (pp. 21)
    std::vector<BigInt> vec_e_power = GenBigIntPowerVector(n,pp.e);

    std::vector<BigInt> vec_zz_temp(LEN);
    for(auto i=0;i<n;i++){
        for(auto j=0;j<pp.RANGE_LEN;j++){
            vec_zz_temp[i*pp.RANGE_LEN+j] = (vec_e_power[i] * vec_short_2_power[j]) * z_square % order;
        }
    }    
    BigInt bn_temp1 = BigIntVectorModInnerProduct(vec_1_power, vec_y_power, BigInt(order)); 
    BigInt bn_temp2 = BigIntVectorModInnerProduct(vec_1_power, vec_zz_temp, BigInt(order)); 

    BigInt bn_c0 = z.ModSub(z_square, order); // z-z^2
    bn_temp1 = bn_c0 * bn_temp1 % order; 
    bn_temp2 = z * bn_temp2 % order; 
  
    BigInt delta_yz = bn_temp1.ModSub(bn_temp2, order);  //Eq (39) also see page 21


    // check Eq (72)  
    ECPoint LEFT = pp.g * proof.taux + pp.h * proof.tx;  // LEFT = g^{\taux} h^\hat{t}

    // the intermediate variables used to compute the right value
    std::vector<ECPoint> vec_A; 
    std::vector<BigInt> vec_a;

    // ECPoint RIGHT = ECPointVectorMul(vec_A, vec_a);  // RIGHT = V^{z^2} h^{\delta_yz} T_1^x T_2^{x^2} 
    ECPoint RIGHT = instance.C[0] * z_square + pp.h * delta_yz + proof.T1 * x + proof.T2 * x_square;

    V1 = (LEFT == RIGHT); 
    #ifdef DEBUG
        std::cout << std::boolalpha << "Condition 1 (Aggregating Log Size BulletProof) = " << V1 << std::endl; 
    #endif

    //check Eq (66,67,68) using Inner Product Argument
    InnerProduct::PP ip_pp = InnerProduct::Setup(LEN, false); 

    ip_pp.vec_g.resize(LEN); 
    std::copy(pp.vec_g.begin(), pp.vec_g.begin()+LEN, ip_pp.vec_g.begin()); // ip_pp.vec_g = pp.vec_g

    ip_pp.vec_h.resize(LEN); 
    std::copy(pp.vec_h.begin(), pp.vec_h.begin()+LEN, ip_pp.vec_h.begin()); 
    std::vector<BigInt> vec_y_inverse_power = GenBigIntPowerVector(LEN, y_inverse); // y^nm
    ip_pp.vec_h = ECPointVectorProduct(ip_pp.vec_h, vec_y_inverse_power);  // ip_pp.vec_h = vec_h_new  

    //InnerProduct_Proof ip_proof = proof.ip_proof;
    InnerProduct::Instance ip_instance;
    ip_pp.u = pp.u * e; // u = u^e 
    
    vec_A.resize(2*ip_pp.VECTOR_LEN+4); 
    std::copy(ip_pp.vec_g.begin(), ip_pp.vec_g.end(), vec_A.begin()); 
    std::copy(ip_pp.vec_h.begin(), ip_pp.vec_h.end(), vec_A.begin()+ip_pp.VECTOR_LEN);

    vec_A[2*ip_pp.VECTOR_LEN] = proof.A; 
    vec_A[2*ip_pp.VECTOR_LEN+1] = proof.S; 
    vec_A[2*ip_pp.VECTOR_LEN+2] = pp.h; 
    vec_A[2*ip_pp.VECTOR_LEN+3] = ip_pp.u; 

    vec_a.resize(2*ip_pp.VECTOR_LEN+4);
    
    std::vector<BigInt> vec_z_minus_unary(LEN, z_minus); 
    std::move(vec_z_minus_unary.begin(), vec_z_minus_unary.end(), vec_a.begin()); // LEFT += g^{-1 z^n} 

    std::vector<BigInt> vec_rr = BigIntVectorModScalar(vec_y_power, z, BigInt(order)); // z y^nm

    for(auto i=0;i<LEN;i++){
        vec_rr[i] = (vec_rr[i] + vec_zz_temp[i]) % order;
    }
    std::move(vec_rr.begin(), vec_rr.end(), vec_a.begin()+ip_pp.VECTOR_LEN); 
     
    vec_a[2*ip_pp.VECTOR_LEN] = bn_1; 
    vec_a[2*ip_pp.VECTOR_LEN+1] = x; 
    vec_a[2*ip_pp.VECTOR_LEN+2] = -proof.mu; 
    vec_a[2*ip_pp.VECTOR_LEN+3] = proof.tx; 


    ip_instance.P = ECPointVectorMul(vec_A, vec_a);  // set P_new = A + S^x + h^{-mu} u^tx  

    V2 = InnerProduct::FastVerify(ip_pp, ip_instance, transcript_str, proof.ip_proof); 
    #ifdef DEBUG
        std::cout << std::boolalpha << "Condition 2 (Aggregating Log Size BulletProof) = " << V2 << std::endl; 
    #endif

    Validity = V1 && V2;     
    #ifdef DEBUG
    if (Validity){ 
        std::cout<< "log size BulletProof accepts >>>" << std::endl; 
    }
    else{
        std::cout<< "log size BulletProof rejects >>>" << std::endl; 
    }
    #endif

    return Validity; 
}

}
#endif
