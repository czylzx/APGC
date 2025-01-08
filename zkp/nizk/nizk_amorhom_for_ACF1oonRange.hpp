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
        ECPoint V; 
        std::vector<ECPoint> vec_c;
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
        fout << proof.eqmdl_proof;   
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
        transcript_str += instance.V.ToByteString();
        for(auto i=0;i<pp.n;i++){
            transcript_str += instance.vec_c[i].ToByteString();
        }
        BigInt p = Hash::StringToBigInt(transcript_str);

        ECPoint init;
        init.SetInfinity();

        size_t n = pp.n;
        std::vector<BigInt> vec_p = GenBigIntPowerVector(2*n, p);
        std::vector<ECPoint> vec_base_new(n*4,init);
        for(auto j=0;j<n;j++){
            vec_base_new[j] = instance.vec_c[j].Invert() * vec_p[2*j] + instance.V.Invert() * vec_p[2*j+1];
            vec_base_new[j+3*n] = pp.g * vec_p[2*j+1];
        }
        
        size_t index = 0;
        for(auto j=0;j<n;j++)
        {
            vec_base_new[n+index] = pp.g * vec_p[2*j];
            vec_base_new[n+index+1] = pp.h * ((vec_p[2*j] + vec_p[2*j+1]) % order); 
            index = index + 2;
        }

        ECPoint F_prime;
        F_prime = ECPointVectorMul(vec_base_new, witness.vec_y);

        // pick random
        std::vector<BigInt> vec_a_init = GenRandomBigIntVectorLessThan(n*4, order);
        BigInt b_inint = GenRandomBigIntLessThan(order);

        // compute Ap,Af
        proof.Ap = ECPointVectorMul(pp.vec_g, vec_a_init) + pp.u * b_inint;
        proof.Af = ECPointVectorMul(vec_base_new, vec_a_init);

        transcript_str = "";
        transcript_str += proof.Ap.ToByteString() + proof.Af.ToByteString();

        //compute the challenge e
        BigInt e = Hash::StringToBigInt(transcript_str);
        
        // compute vec_z
        std::vector<BigInt> vec_z_tmp = BigIntVectorModScalar(witness.vec_y, e, order);
        std::vector<BigInt> vec_z = BigIntVectorModAdd(vec_z_tmp, vec_a_init, order);

        // compute f
        proof.f = (b_inint + e * witness.rp) % order;

        // run eqmdl
        EqmdlProduct2:: PP eqmdl_pp = EqmdlProduct2::Setup(n*4, false);
        eqmdl_pp.vec_g = vec_base_new;
        eqmdl_pp.vec_p = pp.vec_g;

        EqmdlProduct2::Instance eqmdl_instance;
        eqmdl_instance.P = instance.P * e + proof.Ap - (pp.u * proof.f);
        eqmdl_instance.G = F_prime * e + proof.Af;

        EqmdlProduct2::Witness eqmdl_witness;
        eqmdl_witness.vec_a = vec_z;

        transcript_str = "";
        EqmdlProduct2::Prove(eqmdl_pp, eqmdl_instance, eqmdl_witness, transcript_str, proof.eqmdl_proof);

        #ifdef DEBUG
        std::cout << "Amorhom Proof Generation Finishes >>>" << std::endl;
        #endif 
    }

    bool Verify(PP &pp, Instance &instance,  std::string &transcript_str,Proof &proof)
    {
        transcript_str = "";
        transcript_str += instance.P.ToByteString();
        transcript_str += instance.V.ToByteString();
        for(auto i=0;i<pp.n;i++){
            transcript_str += instance.vec_c[i].ToByteString();
        }
        BigInt p = Hash::StringToBigInt(transcript_str);

        ECPoint init;
        init.SetInfinity();

        size_t n = pp.n;
        std::vector<BigInt> vec_p = GenBigIntPowerVector(2*n, p);
        std::vector<ECPoint> vec_base_new(n*4,init);
        for(auto j = 0; j < n; j++)
        {
            vec_base_new[j] = instance.vec_c[j].Invert() * vec_p[2*j] + instance.V.Invert() * vec_p[2*j+1];
            vec_base_new[j+3*n] = pp.g * vec_p[2*j+1];
        }
        
        size_t index = 0;
        for(auto j=0;j<n;j++)
        {
            vec_base_new[n+index] = pp.g * vec_p[2*j];
            vec_base_new[n+index+1] = pp.h * ((vec_p[2*j] + vec_p[2*j+1]) % order); 
            index = index + 2;
        }

        transcript_str = "";
        transcript_str += proof.Ap.ToByteString() + proof.Af.ToByteString();

        //compute the challenge e
        BigInt e = Hash::StringToBigInt(transcript_str);

        EqmdlProduct2:: PP eqmdl_pp = EqmdlProduct2::Setup(n*4, false);
        eqmdl_pp.vec_g = vec_base_new;
        eqmdl_pp.vec_p = pp.vec_g;

        EqmdlProduct2::Instance eqmdl_instance;
        eqmdl_instance.P = instance.P * e + proof.Ap - (pp.u * proof.f);
        eqmdl_instance.G = proof.Af;

        transcript_str = "";
        bool Validity = EqmdlProduct2::Verify(eqmdl_pp, eqmdl_instance, transcript_str, proof.eqmdl_proof);

        #ifdef DEBUG
        if (Validity){ 
            std::cout << "NIZK proof for Amorhom accepts >>>" << std::endl; 
        } else {
            std::cout << "NIZK proof for Amorhom rejects >>>" << std::endl; 
        }
        #endif

        return Validity;

    } 
}
#endif