
#ifndef SOLVENT_EQUAL_HPP
#define SOLVENT_EQUAL_HPP

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"
#include "../nizk/nizk_lin_bit.hpp"
#include "./apgc_amorhom_proof.hpp"
#include "../nizk/nizk_kbit.hpp"
#include "../nizk/nizk_zkdl.hpp"

namespace Solvent_Equal
{

    using Serialization::operator<<;
    using Serialization::operator>>;

    // define the structure of Solvent_Equal Proof
    struct PP
    {
        size_t VECTOR_LEN;     // denotes the size of witness (witness is upto l = 2^VECTOR_LEN)
        size_t LOG_VECTOR_LEN; 
        // size of the vector = VECTOR_LEN
        ECPoint g;
        ECPoint h;
        ECPoint u;
        ECPoint u_new;
        std::vector<ECPoint> vec_g;
        std::vector<ECPoint> vec_g1;
        std::vector<ECPoint> vec_g2;
        std::vector<ECPoint> vec_g3;
        std::vector<ECPoint> vec_g4;
        std::vector<ECPoint> vec_h;

        // ECPoint u;
    };

    // P = vec_g^vec_a vec_h^vec_b u^<vec_a, vec_b>
    struct Instance
    {   
        ECPoint B;
        ECPoint Refresh_CL;
        ECPoint Refresh_CR;
        std::vector<ECPoint> Sum_CL;
        std::vector<ECPoint> Sum_CR;
        std::vector<ECPoint> pk;
    };

    struct Witness
    {
        BigInt sk;
        size_t l0;
        BigInt rb_l0;
        BigInt r_refresh;
        BigInt balance_sender;

    };

    struct Proof
    {
        ECPoint P_equal;
        ECPoint S_equal;
        std::vector<ECPoint> vec_C;
        BigInt zs;
        AmorHom::Proof amorhom_proof;
        LinBit::Proof linbit_proof;
        Kbit::Proof kbit_proof;
        ZkdlProduct::Proof zkdl_product_proof;
        BigInt f_zkdl;
        ECPoint Ax;
    }; 

    std::ofstream &operator<<(std::ofstream &fout, const Proof &proof)
    {
        fout << proof.P_equal << proof.S_equal << proof.vec_C << proof.zs;
        fout << proof.amorhom_proof << proof.linbit_proof << proof.kbit_proof << proof.zkdl_product_proof << proof.f_zkdl << proof.Ax;
        return fout;
    }
    
    std::ifstream &operator>>(std::ifstream &fin,const Proof &proof)
    {
        fin >> proof.P_equal >> proof.S_equal >> proof.vec_C >> proof.zs;
        fin >> proof.amorhom_proof >>proof.linbit_proof >> proof.kbit_proof >> proof.zkdl_product_proof >> proof.f_zkdl >> proof.Ax;
        return fin;
    }
    
    /* (Protocol 2 on pp.15) */
    PP Setup(ECPoint g, ECPoint h, size_t VECTOR_LEN)
    {
        PP pp;
        pp.VECTOR_LEN = VECTOR_LEN;
        pp.LOG_VECTOR_LEN = log2(VECTOR_LEN);
        pp.g = g;
        pp.h = h;
        pp.u = GenRandomGenerator();
        pp.u_new = GenRandomGenerator();
        pp.vec_g = GenRandomECPointVector(pp.LOG_VECTOR_LEN);
        pp.vec_g1 = GenRandomECPointVector(VECTOR_LEN);
        pp.vec_g2 = GenRandomECPointVector(VECTOR_LEN);
        pp.vec_g3 = GenRandomECPointVector(VECTOR_LEN);
        pp.vec_g4 = GenRandomECPointVector(VECTOR_LEN);
        pp.vec_h = GenRandomECPointVector(VECTOR_LEN);
        return pp;
    }

    size_t GetTheNthBit(size_t index, size_t n)
    {
        return (index>>n)&1;
    }

    std::vector<BigInt> lagrange(std::vector<BigInt> &x, std::vector<BigInt> &y) 
    {
        size_t n = x.size();
        if (n != y.size() || n == 0) {
            std::cerr << "Invalid input vectors!" << std::endl;
            exit(EXIT_FAILURE);
        }

        std::vector<BigInt> result(n, bn_0);

        for (size_t i = 0; i < n; i++) {
            BigInt temp = 1;
            for (size_t j = 0; j < n; j++) {
                if (i != j) {
                    temp *= (x[i] - x[j]) % order; 
                }
            }

            BigInt temp_inv = temp.ModInverse(order);
            BigInt L_i = (y[i] * temp_inv) % order; 

            for (size_t j = 0; j < n; j++) {
                BigInt term = L_i;
                for (size_t k = 0; k < n; k++) {
                    if (k != i) {
                        term = (term * (x[j] - x[k]) % order) % order;
                    }
                }
                result[j] = (result[j] + term) % order; 
            }
        }

        return result;
    }

    void Prove(PP pp, Instance instance, Witness witness, std::string &transcript_str, Proof &proof)
    {
        size_t n = pp.VECTOR_LEN;
        size_t m = pp.LOG_VECTOR_LEN;
        std::vector<BigInt> vec_l0;
        for (auto i = 0; i < m; i++)
        {
            if (GetTheNthBit(witness.l0, i) == 1)
            {
                vec_l0.push_back(bn_1);
            }
            else
            {
                vec_l0.push_back(bn_0);
            }
        }
        //LinBit::PP linbit_pp = LinBit::Setup(pp.vec_g, pp.h, n);
        LinBit::PP linbit_pp = LinBit::Setup(n);
        linbit_pp.vec_g = pp.vec_g;
        linbit_pp.h = pp.h;
        LinBit::Instance linbit_instance;
        linbit_instance.P = instance.B;
        
        LinBit::Witness linbit_witness;
        linbit_witness.r = witness.rb_l0;
        linbit_witness.vec_a = vec_l0;


        LinBit::Proof linbit_proof;
        transcript_str = "";
        linbit_proof = LinBit::Prove(linbit_pp, linbit_instance, linbit_witness, transcript_str);

        proof.linbit_proof = linbit_proof;

        std::vector<BigInt> rho = GenRandomBigIntVectorLessThan(m, order);

        //compute p_j_d
        std::vector<std::vector<BigInt>> P; 
        std::vector<BigInt> vec_alpha = linbit_proof.vec_a;
        for(auto j = 0; j < n; j++)
        {
            std::vector<std::vector<BigInt>> A(m, std::vector<BigInt>(2));
            for(auto d = 0; d < m; d++)
            {
                if(GetTheNthBit(j, d) == 1)
                {
                    A[d][0] = vec_alpha[d];
                    A[d][1] = vec_l0[d];
                }
                else
                {
                    A[d][0] = (bn_0 - vec_alpha[d] + order) % order;
                    A[d][1] = (bn_1 - vec_l0[d] + order) % order;
                }
            }
            std::vector<BigInt> p_j = PolyMul(A);
            P.push_back(p_j);
        }
        
        //compute C_deg

        std::vector<ECPoint> vec_C(m);
        for(auto d = 0; d < m; d++)
        {
            vec_C[d] = (pp.vec_h[0] + pp.vec_g1[0]) * P[0][d];
            for(auto j = 1; j < n; j++)
            {
                vec_C[d] += (pp.vec_h[j] + pp.vec_g1[j]) * P[j][d];
            }
            vec_C[d] += pp.u * rho[d];
        }

        //compute s;
        std::vector<BigInt> vec_s(n);
        for(auto j = 0; j < n; j++)
        {
            if(j == witness.l0)
            {
                vec_s[j] = bn_1;
            }
            else
            {
                vec_s[j] = bn_0;
            }
        }

        std::vector<BigInt> vec_sk(n);

        for(auto i =0 ;i < n; i++)
        {
            vec_sk[i] = vec_s[i] * witness.sk % order;
        }

        std::vector<BigInt> vec_r(n);
        for(auto i = 0; i < n; i++)
        {
            vec_r[i] = vec_s[i] * witness.r_refresh % order;
        }
        std::vector<BigInt> vec_v(n);
        for(auto i = 0; i < n; i++)
        {
            vec_v[i] = vec_s[i] * witness.balance_sender % order;
        }

        std::vector<BigInt> vec_y(4*n);
        
        std::copy(vec_s.begin(), vec_s.end(), vec_y.begin());

        std::copy(vec_sk.begin(), vec_sk.end(), vec_y.begin()+n);
        std::copy(vec_r.begin(), vec_r.end(), vec_y.begin()+2*n);
        std::copy(vec_v.begin(), vec_v.end(), vec_y.begin()+3*n);

        //PrintBigIntVector(vec_y, "vec_y");

        BigInt rP = GenRandomBigIntLessThan(order);
        BigInt rS = GenRandomBigIntLessThan(order);

        std::vector<BigInt> vec_1_power(n, bn_1);
        std::vector<BigInt> vec_s_inverse = BigIntVectorModSub(vec_s, vec_1_power, order);
        ECPoint P_equal = ECPointVectorMul(pp.vec_g1, vec_s) + 
                          ECPointVectorMul(pp.vec_g2, vec_sk) + 
                          ECPointVectorMul(pp.vec_g3, vec_r) +
                          ECPointVectorMul(pp.vec_g4, vec_v) +
                          pp.u * rP;
        ECPoint S_equal = ECPointVectorMul(pp.vec_g1, vec_s) + 
                          ECPointVectorMul(pp.vec_h, vec_s_inverse) + 
                          pp.u * rS;

        transcript_str = "";
        transcript_str = transcript_str + linbit_proof.A.ToByteString() + linbit_proof.C.ToByteString() + linbit_proof.D.ToByteString();

        //compute the  challenge x
        BigInt x = Hash::StringToBigInt(transcript_str);

        std::vector<BigInt> vec_p_eval = GenRandomBigIntVectorLessThan(n,order); 

        for(auto j = 0; j < n; j++)
        {
            vec_p_eval[j] = bn_1;
            //std::vector<size_t> vec_index = Decompose(j, 2, m); 
    
            for(auto b = 0; b < m; b++)
            {    
                if(GetTheNthBit(j, b) == 1){
                    vec_p_eval[j] = vec_p_eval[j] * linbit_proof.vec_f[b] % order;
                }
                else{
                    vec_p_eval[j] = vec_p_eval[j] * ((x - linbit_proof.vec_f[b] + order) % order) % order;
                }   
            } 
        }

        BigInt zs = bn_0;
        BigInt x_pow = bn_1;
        for(auto d = 0; d < m; d++)
        {
            zs += rho[d] * x_pow % order;
            x_pow = x_pow * x % order;
        }
        zs = zs + rS * x_pow % order;
        zs = zs % order;

        proof.P_equal = P_equal;
        proof.S_equal = S_equal;
        proof.vec_C = vec_C;
        proof.zs = zs;
        
        std::vector<BigInt> exp_x(m+1);
        exp_x[0] = bn_1;  
        for(auto k = 1; k <= m; k++){
            exp_x[k] = exp_x[k-1] * x % order; 
        }
        
        ECPoint right = proof.vec_C[0] * (bn_0 - exp_x[0] + order);
        for(auto d = 1;d < m; d++){
            right += proof.vec_C[d] * (bn_0 - exp_x[d] + order);
        }

        Kbit::PP kbit_pp = Kbit::Setup(n);
        kbit_pp.g = pp.g;
        kbit_pp.h = pp.h;
        kbit_pp.u = pp.u;
        kbit_pp.u_new = pp.u_new;
        kbit_pp.vec_h = pp.vec_h;
        kbit_pp.vec_g = pp.vec_g1;

        Kbit::Instance kbit_instance;
        kbit_instance.P = S_equal;
        kbit_instance.k = bn_1;

        Kbit::Witness kbit_witness;
        kbit_witness.vec_a = vec_s;
        kbit_witness.vec_b = vec_s_inverse;
        kbit_witness.r = rS;

        Kbit::Proof kbit_proof;
        transcript_str = "";
        kbit_proof = Kbit::Prove(kbit_pp, kbit_instance, kbit_witness, transcript_str);
        proof.kbit_proof = kbit_proof;

        std::vector<BigInt> vec_x_power_m(n,-exp_x[m]);
        ECPoint S_prime = ECPointVectorMul(pp.vec_h, vec_x_power_m);
        S_prime = S_prime + right;
        for(auto i = 0; i < n; i++)
        {
            S_prime = S_prime + (pp.vec_g1[i] + pp.vec_h[i]) * vec_p_eval[i];
        }
        ZkdlProduct::PP zkdl_pp = ZkdlProduct::Setup(4*n, true); //notice the flag
        std::vector<ECPoint> vec_g4zkdl(4*n);
        std::vector<ECPoint> vec_g4amorhom(4*n);
        std::vector<BigInt> vec_a4zkdl(4*n);
        for(auto i = 0; i < n; i++)
        {
            vec_g4zkdl[i] = pp.vec_g2[i];
            vec_g4zkdl[n+i] = pp.vec_g3[i];
            vec_g4zkdl[2*n+i] = pp.vec_g4[i];
            vec_g4zkdl[3*n+i] = pp.vec_h[i];
            vec_a4zkdl[i] = vec_sk[i];
            vec_a4zkdl[n+i] = vec_r[i];
            vec_a4zkdl[2*n+i] = vec_v[i];
            vec_a4zkdl[3*n+i] = -vec_s_inverse[i];
        }
  
        std::copy(vec_g4zkdl.begin(), vec_g4zkdl.end(), vec_g4amorhom.begin());
        zkdl_pp.vec_g = vec_g4amorhom;

        std::vector<BigInt> a_init = GenRandomBigIntVectorLessThan(4*n, order);
        BigInt b_init = GenRandomBigIntLessThan(order);
        ECPoint Ax = ECPointVectorMul(vec_g4amorhom, a_init) + pp.u * b_init;
        proof.Ax = Ax;
        
        std::vector<BigInt> vec_y4zkdl(4*n);
        proof.f_zkdl = b_init + x * (rP - rS);
        std::vector<BigInt> vec_y4zkdl_temp = BigIntVectorScalar(vec_a4zkdl, x);
        vec_y4zkdl = BigIntVectorModAdd(a_init, vec_y4zkdl_temp, order);
        ZkdlProduct::Instance zkdl_instance;
        zkdl_instance.G = (P_equal + S_equal.Invert()) * x + Ax;

        ZkdlProduct::Witness zkdl_witness;
        zkdl_witness.vec_a = vec_y4zkdl;

        ZkdlProduct::Proof zkdl_proof;
        transcript_str = "";
        ZkdlProduct::Prove(zkdl_pp, zkdl_instance,zkdl_witness, transcript_str, zkdl_proof);
        proof.zkdl_product_proof = zkdl_proof;



        AmorHom::PP amorhom_pp = AmorHom::Setup(n);
        amorhom_pp.g = pp.g;
        amorhom_pp.h = pp.h;
        amorhom_pp.u = pp.u;
        amorhom_pp.vec_g = vec_g4amorhom;
        for(auto i = 0; i < n; i++)
        {
            amorhom_pp.vec_g[i] = pp.vec_g1[i];
            amorhom_pp.vec_g[n+i] = pp.vec_g2[i];
            amorhom_pp.vec_g[2*n+i] = pp.vec_g3[i];
            amorhom_pp.vec_g[3*n+i] = pp.vec_g4[i];
        }

        AmorHom::Instance amorhom_instance;
        amorhom_instance.P = P_equal;
        amorhom_instance.S = S_equal;
        amorhom_instance.vec_pk = instance.pk;
        amorhom_instance.Sum_CL = instance.Sum_CL;
        amorhom_instance.Sum_CR = instance.Sum_CR;
        amorhom_instance.Refresh_CL = instance.Refresh_CL;
        amorhom_instance.Refresh_CR = instance.Refresh_CR;
        amorhom_instance.x = x;
        amorhom_instance.zs = zs;

        AmorHom::Witness amorhom_witness;
        amorhom_witness.vec_y = vec_y;
        amorhom_witness.rp = rP;

        AmorHom::Proof amorhom_proof;
        transcript_str = "";
        AmorHom::Prove(amorhom_pp, amorhom_instance, amorhom_witness, transcript_str, amorhom_proof);

        proof.amorhom_proof = amorhom_proof;

    }

    bool Verify(PP &pp, Instance &instance, std::string &transcript_str, Proof &proof)
    {
        bool Validity = false;
        transcript_str = "";
        size_t n = pp.VECTOR_LEN;
        size_t m = pp.LOG_VECTOR_LEN;


        LinBit::PP linbit_pp = LinBit::Setup(n);
        linbit_pp.vec_g = pp.vec_g;
        linbit_pp.h = pp.h;
        LinBit::Instance linbit_instance;
        linbit_instance.P = instance.B;

        LinBit::Proof linbit_proof;
        linbit_proof = proof.linbit_proof;
        transcript_str = "";
        bool V1 = LinBit::Verify(linbit_pp, linbit_instance, transcript_str, linbit_proof);

        transcript_str = "";
        transcript_str = transcript_str + linbit_proof.A.ToByteString() + linbit_proof.C.ToByteString() + linbit_proof.D.ToByteString();

        //compute the  challenge x
        BigInt x = Hash::StringToBigInt(transcript_str);

        std::vector<BigInt> vec_p_eval = GenRandomBigIntVectorLessThan(n,order); 

        for(auto j = 0; j < n; j++)
        {
            vec_p_eval[j] = bn_1;    
            for(auto b = 0; b < m; b++)
            {    
                if(GetTheNthBit(j, b) == 1){
                    vec_p_eval[j] = vec_p_eval[j] * linbit_proof.vec_f[b] % order;
                }
                else{
                    vec_p_eval[j] = vec_p_eval[j] * ((x - linbit_proof.vec_f[b] + order) % order) % order;
                }  
            } 
        }

        std::vector<BigInt> exp_x(m+1);
        exp_x[0] = bn_1;  
        for(auto k = 1; k <= m; k++){
            exp_x[k] = exp_x[k-1] * x % order; 
        }
        ECPoint right = proof.vec_C[0] * (bn_0 - exp_x[0] + order);
        for(auto d = 1;d < m; d++){
            right += proof.vec_C[d] * (bn_0 - exp_x[d] + order);
        }
        Kbit::PP kbit_pp = Kbit::Setup(n);
        kbit_pp.g = pp.g;
        kbit_pp.h = pp.h;
        kbit_pp.u = pp.u;
        kbit_pp.u_new = pp.u_new;
        kbit_pp.vec_h = pp.vec_h;
        kbit_pp.vec_g = pp.vec_g1;

        Kbit::Instance kbit_instance;
        kbit_instance.P = proof.S_equal;
        kbit_instance.k = bn_1;

        Kbit::Proof kbit_proof = proof.kbit_proof;
        transcript_str = "";
        bool V2 = Kbit::Verify(kbit_pp, kbit_instance, transcript_str, kbit_proof);

        std::vector<BigInt> vec_x_power_m(n,-exp_x[m]);
        ECPoint S_prime = ECPointVectorMul(pp.vec_h, vec_x_power_m);
        S_prime = S_prime + right;
        for(auto i = 0; i < n; i++)
        {
            S_prime = S_prime + (pp.vec_g1[i] + pp.vec_h[i]) * vec_p_eval[i];
        }

        ZkdlProduct::PP zkdl_pp = ZkdlProduct::Setup(4*n, true); //notice the flag
        std::vector<ECPoint> vec_g4zkdl(4*n);
        std::vector<ECPoint> vec_g4amorhom(4*n);
        for(auto i = 0; i < n; i++)
        {
            vec_g4zkdl[i] = pp.vec_g2[i];
            vec_g4zkdl[n+i] = pp.vec_g3[i];
            vec_g4zkdl[2*n+i] = pp.vec_g4[i];
            vec_g4zkdl[3*n+i] = pp.vec_h[i];
        }
        std::copy(vec_g4zkdl.begin(), vec_g4zkdl.end(), vec_g4amorhom.begin());
        zkdl_pp.vec_g = vec_g4zkdl;
        ZkdlProduct::Instance zkdl_instance;
        zkdl_instance.G = (proof.P_equal + proof.S_equal.Invert()) * x + proof.Ax - pp.u * proof.f_zkdl;

        ZkdlProduct::Proof zkdl_proof = proof.zkdl_product_proof;
        transcript_str = "";
        bool V3 = ZkdlProduct::Verify(zkdl_pp, zkdl_instance, transcript_str, zkdl_proof);
       
        AmorHom::PP amorhom_pp = AmorHom::Setup(n);
        amorhom_pp.g = pp.g;
        amorhom_pp.h = pp.h;
        amorhom_pp.u = pp.u;
        amorhom_pp.vec_g = vec_g4amorhom;
        for(auto i = 0; i < n; i++)
        {
            amorhom_pp.vec_g[i] = pp.vec_g1[i];
            amorhom_pp.vec_g[n+i] = pp.vec_g2[i];
            amorhom_pp.vec_g[2*n+i] = pp.vec_g3[i];
            amorhom_pp.vec_g[3*n+i] = pp.vec_g4[i];
        }


        AmorHom::Instance amorhom_instance;
        amorhom_instance.P = proof.P_equal;
        amorhom_instance.S = proof.S_equal;
        amorhom_instance.vec_pk = instance.pk;
        amorhom_instance.Sum_CL = instance.Sum_CL;
        amorhom_instance.Sum_CR = instance.Sum_CR;
        amorhom_instance.Refresh_CL = instance.Refresh_CL;
        amorhom_instance.Refresh_CR = instance.Refresh_CR;
        amorhom_instance.x = x;
        amorhom_instance.zs = proof.zs;

        AmorHom::Proof amorhom_proof;
        amorhom_proof = proof.amorhom_proof;
        transcript_str = "";
        bool V4 = AmorHom::Verify(amorhom_pp, amorhom_instance, transcript_str, amorhom_proof);

        ECPoint LEFT = proof.S_equal * (exp_x[m]) + S_prime.Invert();
        bool V5 = (LEFT == pp.u * proof.zs);

        std::cout << "V1 = " << V1 << std::endl;
        std::cout << "V2 = " << V2 << std::endl;
        std::cout << "V3 = " << V3 << std::endl;
        std::cout << "V4 = " << V4 << std::endl;
        std::cout << "V5 = " << V5 << std::endl;
        Validity = V1 && V2 && V3 && V4 && V5;
        return Validity;
    }

}

#endif
