
#ifndef ACF_KOON_HPP
#define ACF_KOON_HPP

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"
#include "../nizk/nizk_lin_bit.hpp"
#include "../compare/nizk_amorhom_for_ACFkoon.hpp"
#include "../nizk/nizk_kbit.hpp"
#include "../nizk/nizk_zkdl.hpp"

namespace ACF_koon
{

    using Serialization::operator<<;
    using Serialization::operator>>;

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

    // define the structure of Solvent_Equal Proof
    struct PP
    {
        size_t VECTOR_LEN;     // denotes the size of witness (witness is upto l = 2^VECTOR_LEN)

        ECPoint g;
        ECPoint h;
        ECPoint u;
        ECPoint u_new;
        std::vector<ECPoint> vec_g;
        std::vector<ECPoint> vec_h;
    };

    struct Instance
    {   
        size_t k;
        std::vector<ECPoint> vec_c;
    };

    struct Witness
    {
        std::vector<BigInt> vec_s;
        std::vector<BigInt> vec_r;
    };

    struct Proof
    {
        ECPoint P;
        ECPoint S;
        std::vector<ECPoint> vec_C;
        BigInt zs;
        AmorHom::Proof amorhom_proof;
        Kbit::Proof kbit_proof;
        BigInt k;
        ZkdlProduct::Proof zkdl_proof;
        BigInt f_zkdl;
        ECPoint Ax;
    }; 

    /* (Protocol 2 on pp.15) */
    PP Setup(size_t VECTOR_LEN)
    {
        PP pp;

        pp.VECTOR_LEN = VECTOR_LEN;
        pp.g = GenRandomGenerator();
        pp.h = GenRandomGenerator();
        pp.u = GenRandomGenerator();
        pp.u_new = GenRandomGenerator();
        pp.vec_g = GenRandomECPointVector(VECTOR_LEN);
        pp.vec_h = GenRandomECPointVector(VECTOR_LEN);

        return pp;
    }


    Proof Prove(PP pp, Instance instance, Witness witness, std::string &transcript_str)
    {
        Proof proof;

        size_t n = pp.VECTOR_LEN;

        // prepare for compute
        std::vector<ECPoint> A;
        std::vector<BigInt> B;

        // compute vec_s o vec_r
        std::vector<BigInt> vec_sor = BigIntVectorProduct(witness.vec_s, witness.vec_r);

        // compute vec_s - vec_1
        std::vector<BigInt> vec_1_n(n, bn_1);
        std::vector<BigInt> vec_s_1 = BigIntVectorModSub(witness.vec_s, vec_1_n, BigInt(order));

        // set vec_y
        std::vector<BigInt> vec_y(2*n,bn_0);
        std::copy(witness.vec_s.begin(), witness.vec_s.end(), vec_y.begin());
        std::copy(vec_sor.begin(), vec_sor.end(),vec_y.begin()+n);

        // pick rP,rS
        BigInt rP = GenRandomBigIntLessThan(order);
        BigInt rS = GenRandomBigIntLessThan(order);

        // compute P
        A.resize(2*n+1);
        B.resize(2*n+1);

        std::copy(pp.vec_g.begin(), pp.vec_g.end(), A.begin());
        std::copy(pp.vec_h.begin(), pp.vec_h.end(), A.begin()+n);
        A[2*n] = pp.u;

        std::copy(witness.vec_s.begin(), witness.vec_s.end(), B.begin());
        std::copy(vec_sor.begin(), vec_sor.end(), B.begin()+n);
        B[2*n] = rP;

        proof.P = ECPointVectorMul(A, B);

        // compute S
        std::copy(witness.vec_s.begin(), witness.vec_s.end(), B.begin());
        std::copy(vec_s_1.begin(), vec_s_1.end(), B.begin()+n);
        B[2*n] = rS;  

        proof.S = ECPointVectorMul(A, B);

        // run kbit
        Kbit::PP kbit_pp = Kbit::Setup(n);
        kbit_pp.g = pp.g;
        kbit_pp.h = pp.h;
        kbit_pp.u = pp.u;
        kbit_pp.u_new = pp.u_new;
        kbit_pp.vec_h = pp.vec_h;
        kbit_pp.vec_g = pp.vec_g;

        Kbit::Instance kbit_instance;
        kbit_instance.P = proof.S;
        kbit_instance.k = instance.k;
        proof.k = BigInt(instance.k);

        Kbit::Witness kbit_witness;
        kbit_witness.vec_a = witness.vec_s;
        kbit_witness.vec_b = vec_s_1;
        kbit_witness.r = rS;

        transcript_str = "";
        proof.kbit_proof = Kbit::Prove(kbit_pp, kbit_instance, kbit_witness, transcript_str);


        // run zkdls
        ZkdlProduct::PP zkdl_pp = ZkdlProduct::Setup(n, true);
        zkdl_pp.vec_g = pp.vec_h;
        std::vector<BigInt> a_init = GenRandomBigIntVectorLessThan(n, order);
        BigInt b_init = GenRandomBigIntLessThan(order);
        proof.Ax = ECPointVectorMul(zkdl_pp.vec_g, a_init) + pp.u * b_init;

        BigInt e_zkdl = Hash::StringToBigInt(proof.Ax.ToByteString());

        proof.f_zkdl = b_init + e_zkdl * (rP - rS);

        std::vector<BigInt> vec_a_zkdl(n);
        for(auto i=0;i<n;i++){
            vec_a_zkdl[i] = (vec_sor[i] + bn_1 - witness.vec_s[i]) % order;
        }
        std::vector<BigInt> vec_y_zkdl_temp = BigIntVectorScalar(vec_a_zkdl, e_zkdl);
        std::vector<BigInt> vec_y_zkdl = BigIntVectorModAdd(a_init, vec_y_zkdl_temp, order);

        ZkdlProduct::Instance zkdl_instance;
        zkdl_instance.G = (proof.P + proof.S.Invert()) * e_zkdl + proof.Ax;

        ZkdlProduct::Witness zkdl_witness;
        zkdl_witness.vec_a = vec_y_zkdl;

        transcript_str = "";
        ZkdlProduct::Prove(zkdl_pp, zkdl_instance,zkdl_witness, transcript_str, proof.zkdl_proof);



        // run amorhom
        std::vector<ECPoint> vec_g_amorhom(2*n);
        for(auto i = 0; i < n; i++)
        {
            vec_g_amorhom[i] = pp.vec_g[i];
            vec_g_amorhom[n+i] = pp.vec_h[i];
        }

        AmorHom::PP amorhom_pp = AmorHom::Setup(n);
        amorhom_pp.g = pp.g;
        amorhom_pp.h = pp.h;
        amorhom_pp.u = pp.u;
        amorhom_pp.vec_g = vec_g_amorhom;

        AmorHom::Instance amorhom_instance;
        amorhom_instance.P = proof.P;
        amorhom_instance.vec_c = instance.vec_c;
        // amorhom_instance.S = S_equal;
        // amorhom_instance.vec_pk = instance.pk;
        // amorhom_instance.Sum_CL = instance.Sum_CL;
        // amorhom_instance.Sum_CR = instance.Sum_CR;
        // amorhom_instance.Refresh_CL = instance.Refresh_CL;
        // amorhom_instance.Refresh_CR = instance.Refresh_CR;
        // amorhom_instance.x = x;
        // amorhom_instance.zs = zs;

        AmorHom::Witness amorhom_witness;
        amorhom_witness.vec_y = vec_y;
        amorhom_witness.rp = rP;

        transcript_str = "";
        AmorHom::Prove(amorhom_pp, amorhom_instance, amorhom_witness, transcript_str, proof.amorhom_proof);

        return proof;
    }

    bool Verify(PP &pp, Instance &instance, std::string &transcript_str, Proof &proof)
    {
        bool Validity = false;

        size_t n = pp.VECTOR_LEN;

        // run kbit
        Kbit::PP kbit_pp = Kbit::Setup(n);
        kbit_pp.g = pp.g;
        kbit_pp.h = pp.h;
        kbit_pp.u = pp.u;
        kbit_pp.u_new = pp.u_new;
        kbit_pp.vec_h = pp.vec_h;
        kbit_pp.vec_g = pp.vec_g;

        Kbit::Instance kbit_instance;
        kbit_instance.P = proof.S;
        kbit_instance.k = proof.k;

        transcript_str = "";
        bool V1 = Kbit::Verify(kbit_pp, kbit_instance, transcript_str, proof.kbit_proof);
        std::cout<<"here"<<std::endl;
        // run zkdl
        ZkdlProduct::PP zkdl_pp = ZkdlProduct::Setup(n, true);
        zkdl_pp.vec_g = pp.vec_h;

        ZkdlProduct::Instance zkdl_instance;
        BigInt e_zkdl = Hash::StringToBigInt(proof.Ax.ToByteString());
        zkdl_instance.G = (proof.P + proof.S.Invert()) * e_zkdl + proof.Ax - pp.u * proof.f_zkdl;

        transcript_str = "";
        bool V2 = ZkdlProduct::Verify(zkdl_pp, zkdl_instance, transcript_str, proof.zkdl_proof);

        // run amorhom
        std::vector<ECPoint> vec_g_amorhom(2*n);
        for(auto i = 0; i < n; i++)
        {
            vec_g_amorhom[i] = pp.vec_g[i];
            vec_g_amorhom[n+i] = pp.vec_h[i];
        }

        AmorHom::PP amorhom_pp = AmorHom::Setup(n);
        amorhom_pp.g = pp.g;
        amorhom_pp.h = pp.h;
        amorhom_pp.u = pp.u;
        amorhom_pp.vec_g = vec_g_amorhom;


        AmorHom::Instance amorhom_instance;
        amorhom_instance.P = proof.P;
        amorhom_instance.vec_c = instance.vec_c;
        // amorhom_instance.S = proof.S_equal;
        // amorhom_instance.vec_pk = instance.pk;
        // amorhom_instance.Sum_CL = instance.Sum_CL;
        // amorhom_instance.Sum_CR = instance.Sum_CR;
        // amorhom_instance.Refresh_CL = instance.Refresh_CL;
        // amorhom_instance.Refresh_CR = instance.Refresh_CR;
        // amorhom_instance.x = x;
        // amorhom_instance.zs = proof.zs;

        transcript_str = "";
        bool V3 = AmorHom::Verify(amorhom_pp, amorhom_instance, transcript_str, proof.amorhom_proof);

        std::cout << "V1 = " << V1 << std::endl;
        std::cout << "V2 = " << V2 << std::endl;
        std::cout << "V3 = " << V3 << std::endl;
        Validity = V1 && V2 && V3;
        return Validity;
    }

}

#endif