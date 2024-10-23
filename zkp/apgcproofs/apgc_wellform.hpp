/***********************************************************************************
this hpp implements the inner product proof system
***********************************************************************************/
#ifndef WELLFORM_PROOF_HPP
#define WELLFORM_PROOF_HPP

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"
// #include <twisted_exponential_elgamal.hpp>
// #include <pedersen.hpp>
#include "../../zkp/apgcproofs/apgc_eqmdl_proof.hpp"

namespace WellFormProduct
{

    using Serialization::operator<<;
    using Serialization::operator>>;

    // define the structure of InnerProduct Proof
    struct PP
    {
        size_t VECTOR_LEN; // denotes the size of witness (witness is upto l = 2^VECTOR_LEN)
        std::vector<ECPoint> vec_g;
        ECPoint g, h, u;
    };

    // P = vec_g^vec_a vec_h^vec_b u^<vec_a, vec_b>
    struct Instance
    {
        std::vector<ECPoint> vec_pk;
        std::vector<ECPoint> vec_CL;
        std::vector<ECPoint> vec_CR;
    };

    struct Witness
    {
        // size of the vector = VECTOR_LEN
        std::vector<BigInt> vec_v;
        std::vector<BigInt> vec_r;
    };

    struct Proof
    {
        // size of the vector = LOG_VECTOR_LEN
        ECPoint G;
        // BigInt c;
        ECPoint AG, ACL, ACR;
        // BigInt e;
        // std::vector<BigInt> z;
        BigInt fg, fcr;
        EqmdlProduct::Proof eqmdl_proof;
    };

    std::ofstream &operator<<(std::ofstream &fout, const WellFormProduct::Proof &proof)
    {
        // fout << proof.G << proof.c << proof.AG << proof.ACL << proof.ACR << proof.e << proof.z << proof.fg << proof.fcr;
        return fout;
    }

    std::ifstream &operator>>(std::ifstream &fin, WellFormProduct::Proof &proof)
    {
        // fin >> proof.G >> proof.c >> proof.AG >> proof.ACL >> proof.ACR >> proof.e >> proof.z >> proof.fg >> proof.fcr;
        return fin;
    }

    std::string ProofToByteString(Proof &proof)
    {
        std::string str;
        // str += proof.G.ToByteString() + proof.c.ToByteString() + proof.AG.ToByteString() + proof.ACL.ToByteString() + proof.ACR.ToByteString() + proof.e.ToByteString();
        // for (auto i = 0; i < proof.z.size(); i++)
        // {
        //     str += proof.z[i].ToByteString();
        // }

        str += proof.fg.ToByteString() + proof.fcr.ToByteString();
        return str;
    }

    void PrintWitness(Witness &witness)
    {
        PrintBigIntVector(witness.vec_v, "vec_v");
        PrintBigIntVector(witness.vec_r, "witness.vec_r");
        // PrintBigIntVector(witness.vec_b, "b");
    }

    void PrintInstance(Instance &instance)
    {
        PrintECPointVector(instance.vec_CL, "vec_CL");
        PrintECPointVector(instance.vec_CR, "vec_CR");
        PrintECPointVector(instance.vec_pk, "vec_pk");
    }

    void PrintProof(Proof &proof)
    {
        // PrintBigIntVector(proof.z, "proof.z");
        // proof.c.Print("proof.c");
        // proof.e.Print("proof.e");
        proof.G.Print("proof.G");
        proof.AG.Print("proof.AG");
        proof.ACL.Print("proof.ACL");
        proof.ACR.Print("proof.ACR");
        proof.fg.Print("proof.fg");
        proof.fcr.Print("proof.fcr");
    };

    /* compute the jth bit of a small integer num \in [0, 2^{l-1}] (count from big endian to little endian) */
    uint64_t GetTheNthBitofInt(uint64_t i, size_t j, size_t LEN)
    {
        uint64_t cursor = 1 << (LEN - j - 1); // set cursor = 2^{m-1} = 1||0...0---(m-1)
        if ((i & cursor) != 0)
            return 1;
        else
            return 0;
    }

    /*  Setup algorithm */
    PP Setup(size_t VECTOR_LEN, TwistedExponentialElGamal::PP pp_enc)
    {
        PP pp;
        pp.g = pp_enc.g;
        pp.h = pp_enc.h;
        pp.u = GenRandomGenerator();
        pp.vec_g = GenRandomECPointVector(VECTOR_LEN);

        return pp;
    }

    // generate NIZK proof for Ci = Enc(pki, v; r) i={1,2,3} the witness is (r, v)
    Proof Prove(PP &pp, Instance &instance, Witness &witness, std::string &transcript_str)
    {
        Proof proof;
        // initialize the transcript with instance
        size_t n = pp.VECTOR_LEN;

        // a
        transcript_str += pp.g.ToByteString() + pp.h.ToByteString() + pp.u.ToByteString();
        BigInt a = Hash::StringToBigInt(transcript_str);

        // G = g * r + a * u
        ECPoint G = ECPointVectorMul(pp.vec_g, witness.vec_r);
        G += pp.u * a;

        proof.G = G;

        // c
        BigInt c = Hash::StringToBigInt(G.ToByteString());


        // CL CR

        std::vector<BigInt> vec_c;
        for (auto i = 0; i < n; i++)
        {
            vec_c[i] = c.ModExp(i, BigInt(order));
        }
        ECPoint CL = ECPointVectorMul(instance.vec_CL, vec_c);
        ECPoint CR = ECPointVectorMul(instance.vec_CR, vec_c);

        // AG
        std::vector<BigInt> vec_a(n);
        BigInt b, p;
        ECPoint AG = ECPointVectorMul(pp.vec_g, vec_a) + pp.u * b;

        proof.AG = AG;

        // ACL
        std::vector<ECPoint> vec_pkc(n);
        for (auto i = 0; i < n; i++)
        {
            vec_pkc[i] = instance.vec_pk[i] * vec_c[i];
        }
        ECPoint ACL = ECPointVectorMul(vec_pkc, vec_a);

        proof.ACL = ACL;

        // ACR
        std::vector<ECPoint> vec_gc(n);
        for (auto i = 0; i < n; i++)
        {
            vec_gc[i] = pp.g * vec_c[i];
        }
        ECPoint ACR = ECPointVectorMul(vec_gc, vec_a) + (pp.h * p);

        proof.ACR = ACR;

        // e
        BigInt e = Hash::StringToBigInt(AG.ToByteString() + ACL.ToByteString() + ACR.ToByteString());

        // proof.e = e;

        // vec_z
        std::vector<BigInt> vec_z;
        std::vector<BigInt> vec_re;
        vec_re = BigIntVectorModScalar(witness.vec_r, e, BigInt(order));
        vec_z = BigIntVectorModAdd(vec_a, vec_re, BigInt(order));

        // proof.z = vec_z;

        // fg  fcr
        BigInt fg, fcr, temp;
        fg = b + a * e;
        for (auto i = 0; i < n; i++)
        {
            temp += vec_c[i] * witness.vec_r[i];
        }
        fcr = p + (temp * e);

        proof.fg = fg % order;
        proof.fcr = fcr % order;

        // sub
        ECPoint g;
        ECPoint pkc;
        ECPoint gc;
        for (auto i = 0; i < n; i++)
        {
            g += pp.vec_g[i];
        }

        pkc = ECPointVectorMul(instance.vec_pk,vec_c);
        gc = ECPointVectorMul({pp.g},vec_c);

        EqmdlProduct::PP eqmdl_pp;
        eqmdl_pp.vec_g = {g};
        eqmdl_pp.vec_h = {pkc};
        eqmdl_pp.vec_p = {gc};
        
        

        EqmdlProduct::Instance eqmdl_ins;
        BigInt fg_inverse = proof.fg.ModInverse(order);
        eqmdl_ins.G = proof.AG + proof.G*e + pp.u*fg_inverse; 

        eqmdl_ins.H = proof.ACL + CL*e;
        
        BigInt fcr_inverse = proof.fcr.ModInverse(order);
        eqmdl_ins.P = proof.ACR + CR*e + pp.h*fcr_inverse;

        EqmdlProduct::Witness wit;
        wit.vec_a = vec_a;
        transcript_str = "";
        

        EqmdlProduct::Prove(eqmdl_pp, eqmdl_ins, wit, transcript_str, proof.eqmdl_proof);

#ifdef DEBUG
        PrintProof(proof);
#endif

        return proof;
    }

    // check NIZK proof PI for Ci = Enc(pki, m; r) the witness is (r1, r2, m)
    bool Verify(PP &pp, Instance &instance, std::string &transcript_str, Proof &proof)
    {
        // initialize the transcript with instance
        bool Validity;

        size_t n = pp.VECTOR_LEN;

        // c
        BigInt c = Hash::StringToBigInt(proof.G.ToByteString());

        // CL CR
        std::vector<BigInt> vec_c;
        for (auto i = 0; i < n; i++)
        {
            vec_c[i] = c.ModExp(i, BigInt(order));
        }
        ECPoint CL = ECPointVectorMul(instance.vec_CL, vec_c);
        ECPoint CR = ECPointVectorMul(instance.vec_CR, vec_c);

        // verify
        ECPoint g;
        ECPoint pkc;
        ECPoint gc;
        for (auto i = 0; i < n; i++)
        {
            g += pp.vec_g[i];
        }

        pkc = ECPointVectorMul(instance.vec_pk,vec_c);
        gc = ECPointVectorMul({pp.g},vec_c);

        EqmdlProduct::PP eqmdl_pp;
        eqmdl_pp.vec_g = {g};
        eqmdl_pp.vec_h = {pkc};
        eqmdl_pp.vec_p = {gc};
        
        

        EqmdlProduct::Instance eqmdl_ins;
        // e
        BigInt e = Hash::StringToBigInt(proof.AG.ToByteString() + proof.ACL.ToByteString() + proof.ACR.ToByteString());
        BigInt fg_inverse = proof.fg.ModInverse(order);
        eqmdl_ins.G = proof.AG + proof.G*e + pp.u*fg_inverse; 

        eqmdl_ins.H = proof.ACL + CL*e;
        
        BigInt fcr_inverse = proof.fcr.ModInverse(order);
        eqmdl_ins.P = proof.ACR + CR*e + pp.h*fcr_inverse;

        EqmdlProduct::Proof eqmdl_proof;
        
        // proof = Prove();
        transcript_str = "";
        
        Validity = EqmdlProduct:: Verify(eqmdl_pp, eqmdl_ins, transcript_str, proof.eqmdl_proof);

        return Validity;
    }

}

#endif
