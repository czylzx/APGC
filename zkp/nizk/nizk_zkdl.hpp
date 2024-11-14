/***********************************************************************************
this hpp implements the eqmdl product proof system
***********************************************************************************/
#ifndef ZKDLProduct_HPP
#define ZKDLProduct_HPP

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"

namespace ZkdlProduct
{

    using Serialization::operator<<;
    using Serialization::operator>>;

    // define the structure of ZkdlProduct Proof
    struct PP
    {
        size_t VECTOR_LEN;     // denotes the size of witness (witness is upto l = 2^VECTOR_LEN)
        size_t LOG_VECTOR_LEN; // LOG_VECTOR_LEN = log(VECTOR_LEN)

        // size of the vector = VECTOR_LEN
        std::vector<ECPoint> vec_g;

        // ECPoint u;
    };

    // P = vec_g^vec_a vec_h^vec_b u^<vec_a, vec_b>
    struct Instance
    {
        ECPoint G;
    };

    struct Witness
    {
        // size of the vector = VECTOR_LEN
        std::vector<BigInt> vec_a;
        // std::vector<BigInt> vec_b;
    };

    struct Proof
    {
        // size of the vector = LOG_VECTOR_LEN

        std::vector<ECPoint> vec_GL;
        std::vector<ECPoint> vec_GR;
        std::vector<BigInt> a;
        // BigInt b;
    };

    std::ofstream &operator<<(std::ofstream &fout, const Proof &proof)
    {
        fout << proof.vec_GL << proof.vec_GR ;
        fout << proof.a;
        return fout;
    }

    std::ifstream &operator>>(std::ifstream &fin, Proof &proof)
    {
        fin >> proof.vec_GL >> proof.vec_GR ;
        fin >> proof.a;
        return fin;
    }

    std::string ProofToByteString(Proof &proof)
    {
        std::string str;

        for (auto i = 0; i < proof.vec_GL.size(); i++)
        {
            str += proof.vec_GL[i].ToByteString();
        }
        for (auto i = 0; i < proof.vec_GR.size(); i++)
        {
            str += proof.vec_GR[i].ToByteString();
        }
        for (auto i = 0; i < proof.a.size(); i++)
        {
            str += proof.a[i].ToByteString();
        }
        return str;
    }

    void PrintPP(PP &pp)
    {
        std::cout << "vector length = " << pp.VECTOR_LEN << std::endl;
        std::cout << "log vector length = " << pp.LOG_VECTOR_LEN << std::endl;

        // size of the vector = VECTOR_LEN
        PrintECPointVector(pp.vec_g, "g");
        // pp.u.Print("u");
    }

    void PrintWitness(Witness &witness)
    {
        PrintBigIntVector(witness.vec_a, "a");
        // PrintBigIntVector(witness.vec_b, "b");
    }

    void PrintInstance(Instance &instance)
    {
        instance.G.Print("ip_instance.G");
    }

    void PrintProof(Proof &proof)
    {
        PrintECPointVector(proof.vec_GL, "L");
        PrintECPointVector(proof.vec_GR, "R");
        PrintBigIntVector(proof.a, "a");
        // proof.b.Print("proof.b");
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

    /* assign left or right part of a Zn vector */
    void AssignBigIntVector(std::vector<BigInt> &result, std::vector<BigInt> &vec_a, std::string selector)
    {
        size_t LEN = vec_a.size() / 2;
        std::vector<BigInt>::iterator start_index;
        if (selector == "left")
            start_index = vec_a.begin();
        if (selector == "right")
            start_index = vec_a.begin() + LEN;

        result.assign(start_index, start_index + LEN);
    }

    // assign left or right part of an ECn vector
    void AssignECPointVector(std::vector<ECPoint> &result, std::vector<ECPoint> &vec_g, std::string selector)
    {
        size_t LEN = vec_g.size() / 2;
        std::vector<ECPoint>::iterator start_index;
        if (selector == "left")
            start_index = vec_g.begin();
        if (selector == "right")
            start_index = vec_g.begin() + LEN;

        result.assign(start_index, start_index + LEN);
    }

    /* this module is used to enable fast verification (cf pp.15) */
    void ComputeVectorSS(std::vector<BigInt> &vec_s, std::vector<BigInt> &vec_x, std::vector<BigInt> &vec_x_inverse)
    {
        size_t m = vec_x.size();
        size_t n = vec_s.size(); // int n = pow(2, m);

        // compute s[0], ..., s[i-1]
        // vector<BIGNUM *> vec_s(n);
        uint64_t flag;
        for (auto i = 0; i < n; i++)
        {
            vec_s[i] = BigInt(bn_1); // set s[i] = 1
            for (auto j = 0; j < m; j++)
            {
                flag = GetTheNthBitofInt(i, j, m);
                if (flag == 1)
                {
                    vec_s[i] = (vec_s[i] * vec_x[j]) % order;
                }
                else
                {
                    vec_s[i] = (vec_s[i] * vec_x_inverse[j]) % order;
                }
            }
        }
    }

    /* (Protocol 2 on pp.15) */
    PP Setup(size_t VECTOR_LEN, bool INITIAL_FLAG)
    {
        // std::cout << "1111 " ;
        PP pp;
        if (IsPowerOfTwo(VECTOR_LEN) == false)
        {
            std::cerr << "VECTOR_LEN must be power of 2" << std::endl;
            exit(EXIT_FAILURE);
        }

        pp.VECTOR_LEN = VECTOR_LEN;
        pp.LOG_VECTOR_LEN = log2(VECTOR_LEN);

        // pp.vec_g.resize(VECTOR_LEN);
        // pp.vec_h.resize(VECTOR_LEN);
        // pp.vec_p.resize(VECTOR_LEN);

        if (INITIAL_FLAG == true)
        {
            pp.vec_g = GenRandomECPointVector(pp.VECTOR_LEN);
            // pp.u = GenRandomGenerator();
        }

        return pp;
    }

    /*
        Generate an argument PI for Relation 3 on pp.13: P = g^a h^b u^<a,b>
        transcript_str is introduced to be used as a sub-protocol
    */
    void Prove(PP pp, Instance instance, Witness witness, std::string &transcript_str, Proof &proof)
    {
        // if (pp.vec_g.size() != pp.vec_h.size())
        // {
        //     std::cerr << "vector size does not match!" << std::endl;
        //     exit(EXIT_FAILURE);
        // }

        if (IsPowerOfTwo(pp.VECTOR_LEN) == false)
        {
            std::cerr << "VECTOR_LEN must be power of 2" << std::endl;
            exit(EXIT_FAILURE);
        }

        size_t n = pp.VECTOR_LEN; // the current size of vec_G and vec_H

        // the last round
        if (n == 1)
        {
            proof.a = witness.vec_a;
            // proof.b = witness.vec_b[0];

#ifdef DEBUG
            std::cerr << "Zkdl Product Proof Generation Finishes >>>" << std::endl;
#endif

            return;
        }

        else
        {
            n = n / 2;

            // prepare the log(n)-th round message
            std::vector<BigInt> vec_aL(n), vec_aR(n);
            std::vector<ECPoint> vec_gL(n), vec_gR(n);

            // prepare aL, aR, bL, bR
            AssignBigIntVector(vec_aL, witness.vec_a, "left");
            AssignBigIntVector(vec_aR, witness.vec_a, "right");

            AssignECPointVector(vec_gL, pp.vec_g, "left");
            AssignECPointVector(vec_gR, pp.vec_g, "right");
            // compute GL,GR,HL,HR,PL,PR

            ECPoint GL = ECPointVectorMul(vec_gL, vec_aR);
            ECPoint GR = ECPointVectorMul(vec_gR, vec_aL);


            proof.vec_GL.emplace_back(GL);
            proof.vec_GR.emplace_back(GR);

            // compute the challenge
            transcript_str += GL.ToByteString() + GR.ToByteString() ;
            BigInt x = Hash::StringToBigInt(transcript_str); // compute the n-th round challenge Eq (26,27)

            // x.Print("P:x");
            // generate new pp
            /*
            ** pp_sub.VECTOR_LEN = pp.VECTOR_LEN/2;
            ** pp_sub.LOG_VECTOR_LEN = pp.LOG_VECTOR_LEN - 1;
            */
            // g'=gL^x gR h'=hL^x hR p'=pL^x pR
            PP pp_sub = Setup(pp.VECTOR_LEN / 2, false);

            // std::cout << vec_gL << std::endl;

            // compute vec_g'
            // for (auto i = 0; i < vec_gL.size(); i++)
            // {
            //     std::cout << "i = " << i << std::endl;
            //     vec_gL[i].Print("P:vec_gL");
            // }

            vec_gL = ECPointVectorScalar(vec_gL, x);

            pp_sub.vec_g = ECPointVectorAdd(vec_gL, vec_gR);


            // generate new instance G' H' P'
            Instance instance_sub;

            std::vector<ECPoint> vec_G_A(2 * n + 1);
            std::vector<BigInt> vec_a(2 * n + 1);


            vec_G_A.clear();

            vec_a.clear();

            vec_G_A = {GL, instance.G, GR};
            vec_a = {x.ModSquare(order), x, bn_1};

            // instance_sub.G = GL * x_square+G*x+GR
            instance_sub.G = ECPointVectorMul(vec_G_A, vec_a);

            // generate new witness a'=aL+x*aR
            Witness witness_sub;
            vec_aR = BigIntVectorModScalar(vec_aR, x, BigInt(order));
            witness_sub.vec_a = BigIntVectorModAdd(vec_aL, vec_aR, BigInt(order));

            // ga.Print("V:ga");
            // G.Print("V:G");
            // std::cout << pp_sub.vec_g.size()<< std::endl;

            // for(auto i = 0; i < pp_sub.vec_g.size(); i++)
            // {
            //     pp_sub.vec_g[i].Print("P:veg_g");
            //     std::cout  << std::endl;
            // }

            Prove(pp_sub, instance_sub, witness_sub, transcript_str, proof);
        }
    }

    /* Check if PI is a valid proof for eqmdl product statement (G1^w = H1 and G2^w = H2) */
    bool Verify(PP &pp, Instance &instance, std::string &transcript_str, Proof &proof)
    {
        if (IsPowerOfTwo(pp.VECTOR_LEN) == false)
        {
            std::cerr << "VECTOR_LEN must be power of 2" << std::endl;
            exit(EXIT_FAILURE);
        }

        bool Validity;

        std::vector<BigInt> vec_aL(pp.LOG_VECTOR_LEN), vec_aR(pp.LOG_VECTOR_LEN);
        std::vector<ECPoint> vec_gL(pp.LOG_VECTOR_LEN), vec_gR(pp.LOG_VECTOR_LEN);
        std::vector<ECPoint> vec_G_A(pp.LOG_VECTOR_LEN);


        // prepare aL, aR, bL, bR
        AssignECPointVector(vec_gL, pp.vec_g, "left");
        AssignECPointVector(vec_gR, pp.vec_g, "right");

        ECPoint G = instance.G;

        PP pp_sub;

        for (auto i = 0; i < pp.LOG_VECTOR_LEN; i++)
        {

            // compute the challenge x
            transcript_str += proof.vec_GL[i].ToByteString() + proof.vec_GR[i].ToByteString() ;
            BigInt x = Hash::StringToBigInt(transcript_str);

            // x.Print("V:x");

            // std::cout << "1111 " ;
            // g'=gL^x gR h'=hL^x hR p'=pL^x pR

            // for (auto i = 0; i < vec_gL.size(); ++i)
            // {
            //     std::cout << "  i=" << i << std::endl;
            //     vec_gL[i].Print("V:vec_gL");
            // }

            // compute vec_g'
            vec_gL = ECPointVectorScalar(vec_gL, x);

            pp_sub.vec_g = ECPointVectorAdd(vec_gL, vec_gR);

            // prepare aL, aR, bL, bR
            AssignECPointVector(vec_gL, pp_sub.vec_g, "left");
            AssignECPointVector(vec_gR, pp_sub.vec_g, "right");

            // generate new instance G' H' P'
            Instance instance_sub;

            vec_G_A = {proof.vec_GL[i], G, proof.vec_GR[i]};
            std::vector<BigInt> vec_a(3);
            vec_a = {x.ModSquare(order), x, bn_1};

            // instance_sub.G = GL * x_square+G*x+GR
            G = ECPointVectorMul(vec_G_A, vec_a);
        }

        // std::vector<BigInt> a = {proof.a[0]};
        ECPoint ga = ECPointVectorMul({pp_sub.vec_g[0]}, proof.a);


        if (G == ga )
        {
            Validity = true;
#ifdef DEBUG
           
            std::cout << "Eqmdl Product Proof Accept >>>" << std::endl;
#endif
        }
        else
        {
            Validity = false;
            
#ifdef DEBUG
            
            std::cout << "Eqmdl Product Proof Reject >>>" << std::endl;
#endif
        }

        return Validity;
    }

}


#endif
