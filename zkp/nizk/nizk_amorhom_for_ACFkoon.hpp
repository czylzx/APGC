#ifndef AmorHom_HPP_
#define AmorHom_HPP_

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"
#include "../../commitment/pedersen.hpp"
#include "../../utility/polymul.hpp"
#include "../bulletproofs/innerproduct_proof.hpp"
#include "../apgcproofs/apgc_eqmdl4two_proof.hpp"
#include <vector>
#include <iostream>

namespace AmorHom
{

    using Serialization::operator<<;
    using Serialization::operator>>;

    struct PP
    {
        size_t n;
        size_t m;                          
        ECPoint g;  
        ECPoint h;
        ECPoint u;
        std::vector<ECPoint> vec_g;                 
    };

    struct Instance
    {
        ECPoint P; 
        std::vector<ECPoint> vec_c;
        // ECPoint S;
        // std::vector<ECPoint> vec_pk;
        // std::vector<ECPoint> Sum_CL;
        // std::vector<ECPoint> Sum_CR;
        // ECPoint Refresh_CL;
        // ECPoint Refresh_CR;
        // BigInt x;
        // BigInt zs;
    };

    struct Witness
    {
        std::vector<BigInt> vec_y;
        BigInt rp;
    };

    struct Proof
    {
        EqmdlProduct2::Proof eqmdl_proof;
        ECPoint Ap;
        ECPoint Af;
        BigInt f;
    };


    std::ofstream &operator<<(std::ofstream &fout, const Proof &proof)
    {
        fout << proof.Ap << proof.Af << proof.f; 
        four << proof.eqmdl_proof;   
        return fout;
    }

    std::ifstream &operator>>(std::ifstream &fin, Proof &proof)
    {
        fin >> proof.Ap >> proof.Af >> proof.f; 
        fin >> proof.eqmdl_proof;   
        return fin;
    } 

    PP Setup(size_t n)
    {  
        PP pp;
        pp.n = n;
        pp.m = log2(n);                                 
        return pp;
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

    void Prove(PP &pp, Instance &instance, Witness &witness,  std::string &transcript_str,Proof &proof)
    {
        transcript_str = "";
        transcript_str += instance.P.ToByteString();
        for(auto i = 0; i < pp.n; i++){
            transcript_str += instance.vec_c[i].ToByteString();
        }

        BigInt p = Hash::StringToBigInt(transcript_str);

        ECPoint init;
        init.SetInfinity();

        size_t n = pp.n;
        std::vector<BigInt> vec_p = GenBigIntPowerVector(pp.n*2+1, p);
        std::vector<ECPoint> vec_base_new(n*2,init);
        for(auto j = 0; j < n; j++)
        {
            vec_base_new[j] = instance.vec_c[j] * vec_p[j];
            vec_base_new[j+n] = pp.g.Invert() * vec_p[j];
        }

        ECPoint F_prime;
        F_prime = ECPointVectorMul(vec_base_new, witness.vec_y);
        std::vector<BigInt> vec_a_init = GenRandomBigIntVectorLessThan(pp.n*2, order);
        BigInt b_inint = GenRandomBigIntLessThan(order);

        ECPoint Ap = ECPointVectorMul(pp.vec_g, vec_a_init) + pp.u * b_inint;
        ECPoint Af = ECPointVectorMul(vec_base_new, vec_a_init);

        proof.Af = Af;
        proof.Ap = Ap;

        transcript_str = "";
        transcript_str += Ap.ToByteString() + Af.ToByteString();

        //compute the challenge e
        BigInt e = Hash::StringToBigInt(transcript_str);
        
        std::vector<BigInt> vec_z_tmp = BigIntVectorModScalar(witness.vec_y, e, order);

        std::vector<BigInt> vec_z = BigIntVectorModAdd(vec_z_tmp, vec_a_init, order);

        BigInt f = (b_inint + e * witness.rp) % order;

        proof.f = f; 
        EqmdlProduct2:: PP eqmdl_pp = EqmdlProduct2::Setup(pp.n*2, true);
        eqmdl_pp.vec_g = vec_base_new;
        eqmdl_pp.vec_p = pp.vec_g;

        EqmdlProduct2::Instance eqmdl_instance;
        eqmdl_instance.P = instance.P * e + Ap - (pp.u * f);
        eqmdl_instance.G = F_prime * e + Af;

        EqmdlProduct2::Witness eqmdl_witness;
        eqmdl_witness.vec_a = vec_z;

        EqmdlProduct2::Proof eqmdl_proof;
        transcript_str = "";
        EqmdlProduct2::Prove(eqmdl_pp, eqmdl_instance, eqmdl_witness, transcript_str, eqmdl_proof);

        proof.eqmdl_proof = eqmdl_proof;
    }

    bool Verify(PP &pp, Instance &instance,  std::string &transcript_str,Proof &proof)
    {
        transcript_str = "";
        transcript_str += instance.P.ToByteString();
        for(auto i = 0; i < pp.n; i++){
            transcript_str += instance.vec_c[i].ToByteString();
        }

        BigInt p = Hash::StringToBigInt(transcript_str);
        
        ECPoint init;
        init.SetInfinity();

        size_t n = pp.n;
        std::vector<BigInt> vec_p = GenBigIntPowerVector(pp.n*2+1, p);
        std::vector<ECPoint> vec_base_new(n*2,init);
        for(auto j = 0; j < n; j++)
        {
            vec_base_new[j] = instance.vec_c[j] * vec_p[j];
            vec_base_new[j+n] = pp.g.Invert() * vec_p[j];
        }

        transcript_str = "";
        transcript_str += proof.Ap.ToByteString() + proof.Af.ToByteString();

        //compute the challenge e
        BigInt e = Hash::StringToBigInt(transcript_str);

        EqmdlProduct2:: PP eqmdl_pp = EqmdlProduct2::Setup(pp.n*2, true);
        eqmdl_pp.vec_g = vec_base_new;
        eqmdl_pp.vec_p = pp.vec_g;

        EqmdlProduct2::Instance eqmdl_instance;
        eqmdl_instance.P = instance.P * e + proof.Ap - (pp.u * proof.f);
        eqmdl_instance.G = proof.Af;

        EqmdlProduct2::Proof eqmdl_proof;

        transcript_str = "";
        eqmdl_proof = proof.eqmdl_proof;
        bool Validity = EqmdlProduct2::Verify(eqmdl_pp, eqmdl_instance, transcript_str, eqmdl_proof);

        return Validity;

    } 
}
#endif