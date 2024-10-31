/***********************************************************************************
this hpp implements the inner product proof system
***********************************************************************************/
#ifndef MutiliPlaintextEquality_PROOF_HPP
#define MutiliPlaintextEquality_PROOF_HPP

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"
// #include <twisted_exponential_elgamal.hpp>
// #include <pedersen.hpp>
//#include "../../zkp/apgcproofs/apgc_eqmdl_proof.hpp"

namespace MutiliPlaintextEquality
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
        ECPoint pk_a;
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
        ECPoint X, Y;
        BigInt z, f;
        
    };

    std::ofstream &operator<<(std::ofstream &fout, const MutiliPlaintextEquality::Proof &proof)
    {
        fout << proof.X << proof.Y << proof.z << proof.f;   
        return fout;
    }

    std::ifstream &operator>>(std::ifstream &fin, MutiliPlaintextEquality::Proof &proof)
    {
        fin >> proof.X >> proof.Y >> proof.z >> proof.f;
        return fin;
    }

    std::string ProofToByteString(Proof &proof)
    {
        std::string str = "";
        // str += proof.G.ToByteString() + proof.c.ToByteString() + proof.AG.ToByteString() + proof.ACL.ToByteString() + proof.ACR.ToByteString() + proof.e.ToByteString();
        // for (auto i = 0; i < proof.z.size(); i++)
        // {
        //     str += proof.z[i].ToByteString();
        // }

        str += proof.X.ToByteString() + proof.Y.ToByteString() + proof.z.ToByteString() + proof.f.ToByteString();
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
        instance.pk_a.Print("pk_a");
        PrintECPointVector(instance.vec_CL, "vec_CL");
        PrintECPointVector(instance.vec_CR, "vec_CR");
    }

    void PrintProof(Proof &proof)
    {
        proof.X.Print("X");
        proof.Y.Print("Y");
        proof.z.Print("z");
        proof.f.Print("f");
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
    PP Setup(ECPoint g,ECPoint h, size_t VECTOR_LEN)
    {

        PP pp;
        pp.VECTOR_LEN = VECTOR_LEN;
        pp.g = g;// need to be modified
        pp.h = h;// need to be modified
        //pp.u = GenRandomGenerator();
        //pp.vec_g = GenRandomECPointVector(VECTOR_LEN);

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

    // generate NIZK proof for Ci = Enc(pki, v; r) i={1,2,3} the witness is (r, v)
    void Prove(PP &pp, Instance &instance, Witness &witness, std::string &transcript_str, Proof &proof)
    {
        // Proof proof;
        // initialize the transcript with instance
        size_t n = pp.VECTOR_LEN;

        // rx, mx
        BigInt rx = GenRandomBigIntLessThan(BigInt(order));
        BigInt mx = GenRandomBigIntLessThan(BigInt(order));

        // compute X, Y
        ECPoint X = instance.pk_a * rx;
        ECPoint Y = pp.g * rx + pp.h * mx;

        proof.X = X;
        proof.Y = Y;

        // compute e
        BigInt e = Hash::StringToBigInt(X.ToByteString() + Y.ToByteString());

        // compute z = rx + e * vec_r
        std::vector<BigInt> vec_e_power = GenBigIntPowerVector(n, e);

        BigInt z = rx + BigIntVectorModInnerProduct(witness.vec_r, vec_e_power, BigInt(order));
        BigInt f = mx + BigIntVectorModInnerProduct(witness.vec_v, vec_e_power, BigInt(order));

        proof.z = z;
        proof.f = f;

// #ifdef DEBUG
//         PrintProof(proof);
// #endif

        // return proof;
    }

    bool Verify(PP &pp, Instance &instance, std::string &transcript_str, Proof &proof)
    {
        // initialize the transcript with instance
        bool Validity;

        size_t n = pp.VECTOR_LEN;
        transcript_str = "";

        // compute e
        BigInt e = Hash::StringToBigInt(proof.X.ToByteString() + proof.Y.ToByteString());

        // check
        ECPoint LEFT1 = instance.pk_a * proof.z;
        ECPoint LEFT2 = pp.g * proof.z + pp.h * proof.f;

        bool V1 = false;
        bool V2 = false;
        std::vector<BigInt> vec_e_power = GenBigIntPowerVector(n, e);
        ECPoint RIGHT1 = ECPointVectorMul(instance.vec_CL, vec_e_power) + proof.X;
        ECPoint RIGHT2 = ECPointVectorMul(instance.vec_CR, vec_e_power) + proof.Y;

        if(LEFT1 == RIGHT1)
        {
            V1 = true;
            //std::cout << "V1 is true" << std::endl;
        }
        if(LEFT2 == RIGHT2)
        {
            V2 = true;
            //std::cout << "V2 is true" << std::endl;
        }
        Validity = V1 && V2;
        return Validity;
    }

}

#endif
