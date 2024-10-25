/***********************************************************************************
this hpp implements the inner product proof system
***********************************************************************************/
#ifndef EQ_PROOF_HPP
#define EQ_PROOF_HPP

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"

namespace EqmdlProduct
{

    using Serialization::operator<<;
    using Serialization::operator>>;

    // define the structure of InnerProduct Proof
    struct PP
    {
        size_t VECTOR_LEN;     // denotes the size of witness (witness is upto l = 2^VECTOR_LEN)
        size_t LOG_VECTOR_LEN; // LOG_VECTOR_LEN = log(VECTOR_LEN)

        // size of the vector = VECTOR_LEN
        //std::vector<ECPoint> vec_g{VECTOR_LEN};
        //std::vector<ECPoint> vec_h{VECTOR_LEN};
        //std::vector<ECPoint> vec_p{VECTOR_LEN};
        std::vector<ECPoint> vec_g;
        std::vector<ECPoint> vec_h;
        std::vector<ECPoint> vec_p;
        // ECPoint u;
    };

    // P = vec_g^vec_a vec_h^vec_b u^<vec_a, vec_b>
    struct Instance
    {
        ECPoint P;
        ECPoint G;
        ECPoint H;
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
        std::vector<ECPoint> vec_PL;
        std::vector<ECPoint> vec_PR;
        std::vector<ECPoint> vec_GL;
        std::vector<ECPoint> vec_GR;
        std::vector<ECPoint> vec_HL;
        std::vector<ECPoint> vec_HR;
        std::vector<BigInt> a;
        // BigInt b;
    };

    std::ofstream &operator<<(std::ofstream &fout, const EqmdlProduct::Proof &proof)
    {
        fout << proof.vec_PL << proof.vec_PR << proof.vec_GL << proof.vec_GR << proof.vec_HL << proof.vec_HR;
        fout << proof.a;
        return fout;
    }

    std::ifstream &operator>>(std::ifstream &fin, EqmdlProduct::Proof &proof)
    {
        fin >> proof.vec_PL >> proof.vec_PR >> proof.vec_GL >> proof.vec_GR >> proof.vec_HL >> proof.vec_HR;
        fin >> proof.a;
        return fin;
    }

    std::string ProofToByteString(Proof &proof)
    {
        std::string str;
        for (auto i = 0; i < proof.vec_PL.size(); i++)
        {
            str += proof.vec_PL[i].ToByteString();
        }
        for (auto i = 0; i < proof.vec_PR.size(); i++)
        {
            str += proof.vec_PR[i].ToByteString();
        }

        for (auto i = 0; i < proof.vec_GL.size(); i++)
        {
            str += proof.vec_PL[i].ToByteString();
        }
        for (auto i = 0; i < proof.vec_GR.size(); i++)
        {
            str += proof.vec_PR[i].ToByteString();
        }

        for (auto i = 0; i < proof.vec_HL.size(); i++)
        {
            str += proof.vec_HL[i].ToByteString();
        }
        for (auto i = 0; i < proof.vec_HR.size(); i++)
        {
            str += proof.vec_HR[i].ToByteString();
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
        PrintECPointVector(pp.vec_h, "h");
        PrintECPointVector(pp.vec_p, "p");

        // pp.u.Print("u");
    }

    void PrintWitness(Witness &witness)
    {
        PrintBigIntVector(witness.vec_a, "a");
        // PrintBigIntVector(witness.vec_b, "b");
    }

    void PrintInstance(Instance &instance)
    {
        instance.P.Print("ip_instance.P");
    }

    void PrintProof(Proof &proof)
    {
        PrintECPointVector(proof.vec_PL, "L");
        PrintECPointVector(proof.vec_PR, "R");
        PrintECPointVector(proof.vec_GL, "L");
        PrintECPointVector(proof.vec_GR, "R");
        PrintECPointVector(proof.vec_HL, "L");
        PrintECPointVector(proof.vec_HR, "R");
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
        std::cout << "1111 " ;
        PP pp;
        if (IsPowerOfTwo(VECTOR_LEN) == false)
        {
            std::cerr << "VECTOR_LEN must be power of 2" << std::endl;
            exit(EXIT_FAILURE);
        }

        pp.VECTOR_LEN = VECTOR_LEN;
        pp.LOG_VECTOR_LEN = log2(VECTOR_LEN);

        if (INITIAL_FLAG == true)
        {
            pp.vec_g = GenRandomECPointVector(pp.VECTOR_LEN);
            pp.vec_h = GenRandomECPointVector(pp.VECTOR_LEN);
            pp.vec_p = GenRandomECPointVector(pp.VECTOR_LEN);
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
            std::cerr << "Inner Product Proof Generation Finishes >>>" << std::endl;
#endif

            return;
        }

        else
        {
            n = n / 2;

            // prepare the log(n)-th round message
            std::vector<BigInt> vec_aL(n), vec_aR(n);
            std::vector<ECPoint> vec_gL(n), vec_gR(n), vec_hL(n), vec_hR(n), vec_pL(n), vec_pR(n);

            // prepare aL, aR, bL, bR
            AssignBigIntVector(vec_aL, witness.vec_a, "left");
            AssignBigIntVector(vec_aR, witness.vec_a, "right");
            // AssignBigIntVector(vec_bL, witness.vec_b, "left");
            // AssignBigIntVector(vec_bR, witness.vec_b, "right");

            AssignECPointVector(vec_gL, pp.vec_g, "left");
            AssignECPointVector(vec_gR, pp.vec_g, "right");
            AssignECPointVector(vec_hL, pp.vec_h, "left");
            AssignECPointVector(vec_hR, pp.vec_h, "right");
            AssignECPointVector(vec_pL, pp.vec_p, "left");
            AssignECPointVector(vec_pR, pp.vec_p, "right");

            // compute GL,GR,HL,HR,PL,PR

            ECPoint GL = ECPointVectorMul(vec_gL, vec_aR);
            ECPoint GR = ECPointVectorMul(vec_gR, vec_aL);
            ECPoint HL = ECPointVectorMul(vec_hL, vec_aR);
            ECPoint HR = ECPointVectorMul(vec_hR, vec_aL);
            ECPoint PL = ECPointVectorMul(vec_pL, vec_aR);
            ECPoint PR = ECPointVectorMul(vec_pR, vec_aL);

            proof.vec_GL.emplace_back(GL);
            proof.vec_GR.emplace_back(GR);
            proof.vec_HL.emplace_back(HL);
            proof.vec_HR.emplace_back(HR);
            proof.vec_PL.emplace_back(PL);
            proof.vec_PR.emplace_back(PR);

            // compute the challenge
            transcript_str += GL.ToByteString() + GR.ToByteString() + HL.ToByteString() + HR.ToByteString() + PL.ToByteString() + PR.ToByteString();
            BigInt x = Hash::StringToBigInt(transcript_str); // compute the n-th round challenge Eq (26,27)

            // generate new pp
            /*
            ** pp_sub.VECTOR_LEN = pp.VECTOR_LEN/2;
            ** pp_sub.LOG_VECTOR_LEN = pp.LOG_VECTOR_LEN - 1;
            */
            // g'=gL^x gR h'=hL^x hR p'=pL^x pR
            PP pp_sub = Setup(pp.VECTOR_LEN / 2, false);

            // compute vec_g'
            vec_gL = ECPointVectorScalar(vec_gL, x);
            pp_sub.vec_g = ECPointVectorAdd(vec_gL, vec_gR);

            // compute vec_h'
            vec_hL = ECPointVectorScalar(vec_hL, x);
            pp_sub.vec_h = ECPointVectorAdd(vec_hL, vec_hR);

            // compute vec_p'
            vec_pL = ECPointVectorScalar(vec_pL, x);
            pp_sub.vec_p = ECPointVectorAdd(vec_pL, vec_pR);

            // generate new instance G' H' P'
            Instance instance_sub;

            std::vector<ECPoint> vec_G_A(2 * n + 1);
            std::vector<BigInt> vec_a(2 * n + 1);
            std::vector<ECPoint> vec_H_A(2 * n + 1);
            std::vector<ECPoint> vec_P_A(2 * n + 1);

            vec_G_A = {GL, instance.G, GR};
            vec_H_A = {HL, instance.H, HR};
            vec_P_A = {PL, instance.P, PR};
            vec_a = {x.ModSquare(order), x, bn_1};

            // instance_sub.G = GL * x_square+G*x+GR
            instance_sub.G = ECPointVectorMul(vec_G_A, vec_a);
            instance_sub.H = ECPointVectorMul(vec_H_A, vec_a);
            instance_sub.P = ECPointVectorMul(vec_P_A, vec_a);

            // generate new witness a'=aL+x*aR
            Witness witness_sub;
            vec_aR = BigIntVectorModScalar(vec_aR, x, BigInt(order));
            witness_sub.vec_a = BigIntVectorModAdd(vec_aL, vec_aR, BigInt(order));

            Prove(pp_sub, instance_sub, witness_sub, transcript_str, proof);
        }
    }

    /* Check if PI is a valid proof for inner product statement (G1^w = H1 and G2^w = H2) */
    bool Verify(PP &pp, Instance &instance, std::string &transcript_str, Proof &proof)
    {
        if (IsPowerOfTwo(pp.VECTOR_LEN) == false)
        {
            std::cerr << "VECTOR_LEN must be power of 2" << std::endl;
            exit(EXIT_FAILURE);
        }

        bool Validity;

        std::vector<BigInt> vec_aL(pp.LOG_VECTOR_LEN), vec_aR(pp.LOG_VECTOR_LEN);
        std::vector<ECPoint> vec_gL(pp.LOG_VECTOR_LEN), vec_gR(pp.LOG_VECTOR_LEN), vec_hL(pp.LOG_VECTOR_LEN), vec_hR(pp.LOG_VECTOR_LEN), vec_pL(pp.LOG_VECTOR_LEN), vec_pR(pp.LOG_VECTOR_LEN);
        std::vector<ECPoint> vec_G_A(pp.LOG_VECTOR_LEN );
        std::vector<ECPoint> vec_H_A(pp.LOG_VECTOR_LEN );
        std::vector<ECPoint> vec_P_A(pp.LOG_VECTOR_LEN );

        // prepare aL, aR, bL, bR
        AssignECPointVector(vec_gL, pp.vec_g, "left");
        AssignECPointVector(vec_gR, pp.vec_g, "right");
        AssignECPointVector(vec_hL, pp.vec_h, "left");
        AssignECPointVector(vec_hR, pp.vec_h, "right");
        AssignECPointVector(vec_pL, pp.vec_p, "left");
        AssignECPointVector(vec_pR, pp.vec_p, "right");

        ECPoint G = instance.G, H = instance.H, P = instance.P;

        PP pp_sub;

        for (auto i = 0; i < pp.LOG_VECTOR_LEN; i++)
        {

            // compute the challenge x
            transcript_str += proof.vec_GL[i].ToByteString() + proof.vec_GR[i].ToByteString()
                            + proof.vec_HL[i].ToByteString() + proof.vec_HR[i].ToByteString()  
                            + proof.vec_PL[i].ToByteString() + proof.vec_PR[i].ToByteString();
            BigInt x = Hash::StringToBigInt(transcript_str);

            std::cout << "1111 " ;
            // g'=gL^x gR h'=hL^x hR p'=pL^x pR

            // compute vec_g'
            vec_gL = ECPointVectorScalar(vec_gL, x);
            pp_sub.vec_g = ECPointVectorAdd(vec_gL, vec_gR);

            // compute vec_h'
            vec_hL = ECPointVectorScalar(vec_hL, x);
            pp_sub.vec_h = ECPointVectorAdd(vec_hL, vec_hR);

            // compute vec_p'
            vec_pL = ECPointVectorScalar(vec_pL, x);
            pp_sub.vec_p = ECPointVectorAdd(vec_pL, vec_pR);


            std::cout << "2222 " ;
            // generate new instance G' H' P'
            Instance instance_sub;

            vec_G_A = {proof.vec_GL[i], G, proof.vec_GR[i]};
            vec_H_A = {proof.vec_HL[i], H, proof.vec_HR[i]};
            vec_P_A = {proof.vec_PL[i], P, proof.vec_PR[i]};
            std::vector<BigInt> vec_a(3);
            vec_a = {x.ModSquare(order), x, bn_1};

            // instance_sub.G = GL * x_square+G*x+GR
            G = ECPointVectorMul(vec_G_A, vec_a);
            H = ECPointVectorMul(vec_H_A, vec_a);
            P = ECPointVectorMul(vec_P_A, vec_a);
        }

        ECPoint ga = ECPointVectorMul({pp_sub.vec_g[0]}, {proof.a});
        ECPoint ha = ECPointVectorMul({pp_sub.vec_h[0]}, {proof.a});
        ECPoint pa = ECPointVectorMul({pp_sub.vec_p[0]}, {proof.a});

        if (G == ga && H == ha && P == pa)
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