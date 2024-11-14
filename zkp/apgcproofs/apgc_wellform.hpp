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
    /*struct ConRandom
    {
        BigInt c;
        BigInt e;
    };*/

    std::ofstream &operator<<(std::ofstream &fout, const WellFormProduct::Proof &proof)
    {
        fout << proof.G << proof.AG << proof.ACL << proof.ACR << proof.fg << proof.fcr;
        fout << proof.eqmdl_proof;
        return fout;
    }

    std::ifstream &operator>>(std::ifstream &fin, WellFormProduct::Proof &proof)
    {
        fin >> proof.G >> proof.AG >> proof.ACL >> proof.ACR >> proof.fg >> proof.fcr;
        fin >> proof.eqmdl_proof;
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
        PrintBigIntVector(witness.vec_v, "witness.vec_v");
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
    PP Setup(ECPoint g, ECPoint h, size_t VECTOR_LEN)
    {

        PP pp;
        pp.VECTOR_LEN = VECTOR_LEN;
        pp.g = g;
        pp.h = h;
        pp.u = GenRandomGenerator();
        pp.vec_g = GenRandomECPointVector(VECTOR_LEN);

        return pp;
    }

    // generate NIZK proof for Ci = Enc(pki, v; r) i={1,2,3} the witness is (r, v)
    void Prove(PP &pp, Instance &instance, Witness &witness, std::string &transcript_str, Proof &proof)
    {
        // Proof proof;
        // initialize the transcript with instance
        size_t n = pp.VECTOR_LEN;

        // a
        BigInt a = GenRandomBigIntLessThan(BigInt(order));

        // G = g * r + a * u
        ECPoint G = ECPointVectorMul(pp.vec_g, witness.vec_r) + pp.u * a;

        proof.G = G;

        // c
        transcript_str += G.ToByteString();
        BigInt c = Hash::StringToBigInt(transcript_str);

        // CL CR

        std::vector<BigInt> vec_c(n);
        for (auto i = 0; i < n; i++)
        {
            vec_c[i] = c.ModExp(i, BigInt(order));
        }
        ECPoint CL = ECPointVectorMul(instance.vec_CL, vec_c);
        ECPoint CR = ECPointVectorMul(instance.vec_CR, vec_c);

        // a,b p
        std::vector<BigInt> vec_a = GenRandomBigIntVectorLessThan(n, BigInt(order));
        BigInt b = GenRandomBigIntLessThan(BigInt(order));
        BigInt p = GenRandomBigIntLessThan(BigInt(order));

        // AG
        ECPoint AG = ECPointVectorMul(pp.vec_g, vec_a) + pp.u * b;

        proof.AG = AG;

        // ACL
        std::vector<ECPoint> vec_pkc(n);
        vec_pkc = ECPointVectorProduct(instance.vec_pk, vec_c);
        ECPoint ACL = ECPointVectorMul(vec_pkc, vec_a);

        proof.ACL = ACL;

        // ACR
        std::vector<ECPoint> vec_gc(n);
        // vec_gc = ECPointVectorProduct(pp.g,vec_c);
        for (auto i = 0; i < n; i++)
        {
            vec_gc[i] = pp.g * vec_c[i];
        }
        ECPoint ACR = ECPointVectorMul(vec_gc, vec_a) + (pp.h * p);

        proof.ACR = ACR;

        // e
        transcript_str = "";
        transcript_str += AG.ToByteString() + ACL.ToByteString() + ACR.ToByteString();
        BigInt e = Hash::StringToBigInt(transcript_str);

        // vec_z
        std::vector<BigInt> vec_z(n);
        std::vector<BigInt> vec_re(n);
        vec_re = BigIntVectorModScalar(witness.vec_r, e, BigInt(order));
        vec_z = BigIntVectorModAdd(vec_a, vec_re, BigInt(order));

        // proof.z = vec_z;

        // fg  fcr
        BigInt fg, fcr, temp = bn_0;
        fg = (b + a * e) % BigInt(order);
        for (auto i = 0; i < n; i++)
        {
            temp += vec_c[i] * witness.vec_v[i];
        }
        fcr = (p + (temp * e)) % BigInt(order);

        proof.fg = fg;
        proof.fcr = fcr;

        // sub
        std::vector<ECPoint> eq_vec_g(n);
        std::vector<ECPoint> eq_vec_pkc(n);
        std::vector<ECPoint> eq_vec_gc(n);

        eq_vec_g = pp.vec_g;
        eq_vec_pkc = ECPointVectorProduct(instance.vec_pk, vec_c);
        for (auto i = 0; i < n; i++)
        {
            eq_vec_gc[i] = pp.g * vec_c[i];
        }
        EqmdlProduct::PP eqmdl_pp = EqmdlProduct::Setup(n, true);

        eqmdl_pp.vec_g = eq_vec_g;
        eqmdl_pp.vec_h = eq_vec_pkc;
        eqmdl_pp.vec_p = eq_vec_gc;

        EqmdlProduct::Instance eqmdl_ins;
        ECPoint u_inverse = pp.u.Invert();

        eqmdl_ins.G = AG + G * e + u_inverse * fg;

        eqmdl_ins.H = ACL + CL * e;

        ECPoint h_inverse = pp.h.Invert();
        eqmdl_ins.P = ACR + CR * e + h_inverse * fcr;

        EqmdlProduct::Witness eqmdl_wit;
        eqmdl_wit.vec_a = vec_z;
        transcript_str = "";

        EqmdlProduct::Prove(eqmdl_pp, eqmdl_ins, eqmdl_wit, transcript_str, proof.eqmdl_proof);

#ifdef DEBUG
        PrintProof(proof);
#endif

        // return proof;
    }

    bool Verify(PP &pp, Instance &instance, std::string &transcript_str, Proof &proof)
    {
        // initialize the transcript with instance
        bool Validity;

        size_t n = pp.VECTOR_LEN;

        // c
        BigInt c = Hash::StringToBigInt(proof.G.ToByteString());
        // c.Print("V:c");

        // CL CR
        std::vector<BigInt> vec_c(n);
        // c^j
        for (auto i = 0; i < n; i++)
        {
            vec_c[i] = c.ModExp(i, BigInt(order));
        }
        ECPoint CL = ECPointVectorMul(instance.vec_CL, vec_c);
        ECPoint CR = ECPointVectorMul(instance.vec_CR, vec_c);

        // verify
        std::vector<ECPoint> eq_vec_g(n);
        std::vector<ECPoint> eq_vec_pkc(n);
        std::vector<ECPoint> eq_vec_gc(n);

        eq_vec_g = pp.vec_g;
        eq_vec_pkc = ECPointVectorProduct(instance.vec_pk, vec_c);
        for (auto i = 0; i < n; i++)
        {
            eq_vec_gc[i] = pp.g * vec_c[i];
        }

        EqmdlProduct::PP eqmdl_pp = EqmdlProduct::Setup(n, false);
        // EqmdlProduct::PP eqmdl_pp;
        eqmdl_pp.VECTOR_LEN = n;
        eqmdl_pp.vec_g = eq_vec_g;
        eqmdl_pp.vec_h = eq_vec_pkc;
        eqmdl_pp.vec_p = eq_vec_gc;

        EqmdlProduct::Instance eqmdl_ins;
        // e
        BigInt e = Hash::StringToBigInt(proof.AG.ToByteString() + proof.ACL.ToByteString() + proof.ACR.ToByteString());

        ECPoint u_inverse = pp.u.Invert();
        eqmdl_ins.G = proof.AG + proof.G * e + u_inverse * proof.fg;

        eqmdl_ins.H = proof.ACL + CL * e;

        ECPoint h_inverse = pp.h.Invert();
        eqmdl_ins.P = proof.ACR + CR * e + h_inverse * proof.fcr;

        transcript_str = "";

        Validity = EqmdlProduct::Verify(eqmdl_pp, eqmdl_ins, transcript_str, proof.eqmdl_proof);

        return Validity;
    }

}

#endif
