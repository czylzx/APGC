/***********************************************************************************
this hpp implements the eqmdl product proof system
***********************************************************************************/
#ifndef SOLVENT_EQUAL_HPP
#define SOLVENT_EQUAL_HPP

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"
#include "../nizk/nizk_lin_bit.hpp"
#include "./apgc_amorhom_proof.hpp"

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
        std::vector<ECPoint> vec_g;
        std::vector<ECPoint> vec_h;
        std::vector<ECPoint> vec_u;
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
        //BigInt l0;
        BigInt rb_l0;
        BigInt r_refresh;
        BigInt balance_sender;
        // size of the vector = VECTOR_LEN
        //std::vector<BigInt> vec_a;
        // std::vector<BigInt> vec_b;
    };

    struct Proof
    {
        ECPoint P_equal;
        std::vector<ECPoint> vec_C;
        BigInt zs;
        AmorHom::Proof amorhom_proof;
        LinBit::Proof linbit_proof;
        
    };

    std::ofstream &operator<<(std::ofstream &fout, const Solvent_Equal::Proof &proof)
    {
        
        return fout;
    }

    std::ifstream &operator>>(std::ifstream &fin, Solvent_Equal::Proof &proof)
    {
        
        return fin;
    }

    std::string ProofToByteString(Proof &proof)
    {
        std::string str;
        
        return str;
    }

    void PrintPP(PP &pp)
    {
        std::cout << "vector length = " << pp.VECTOR_LEN << std::endl;
        // size of the vector = VECTOR_LEN
        PrintECPointVector(pp.vec_g, "g");
        PrintECPointVector(pp.vec_h, "h");
        PrintECPointVector(pp.vec_u, "u");

        // pp.u.Print("u");
    }

    void PrintWitness(Witness &witness)
    {
        
        // PrintBigIntVector(witness.vec_b, "b");
    }

    void PrintInstance(Instance &instance)
    {
        //instance.P.Print("ip_instance.P");
    }

    void PrintProof(Proof &proof)
    {
       
    };

 

    /* (Protocol 2 on pp.15) */
    PP Setup(ECPoint g, ECPoint h, size_t VECTOR_LEN)
    {
        PP pp;
        /*if (IsPowerOfTwo(VECTOR_LEN) == false)
        {
            std::cerr << "VECTOR_LEN must be power of 2" << std::endl;
            exit(EXIT_FAILURE);
        }*/

        pp.VECTOR_LEN = VECTOR_LEN;
        pp.LOG_VECTOR_LEN = log2(VECTOR_LEN);
        pp.g = g;
        pp.h = h;
        pp.u = GenRandomGenerator();
        pp.vec_g = GenRandomECPointVector(pp.LOG_VECTOR_LEN);
        pp.vec_h = GenRandomECPointVector(VECTOR_LEN);
        pp.vec_u = GenRandomECPointVector(4*VECTOR_LEN);
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
                temp *= (x[i] - x[j]) % order; // 计算基多项式的分母
            }
        }

        // 模逆运算
        BigInt temp_inv = temp.ModInverse(order);
        BigInt L_i = (y[i] * temp_inv) % order; // L_i(x)

        // 更新结果
        for (size_t j = 0; j < n; j++) {
            BigInt term = L_i;
            for (size_t k = 0; k < n; k++) {
                if (k != i) {
                    term = (term * (x[j] - x[k]) % order) % order; // 计算每个 L_i
                }
            }
            result[j] = (result[j] + term) % order; // 累加到结果
        }
    }

    return result;
}

    // std::vector<BigInt> lagrange(std::vector<BigInt> x, std::vector<BigInt> y)
    // {
    //     size_t n = x.size();
    //     if(n != y.size())
    //     {
    //         std::cerr << "vector size does not match!" << std::endl;
    //         exit(EXIT_FAILURE);
    //     }
    //     std::vector<BigInt> result(n, bn_0);
    //     std::vector<BigInt> Xi(2*n+1,bn_0);
    //     std::vector<BigInt> Fi(2*n+1,bn_0);
    //     std::vector<BigInt> Xi_temp(2*n+1,bn_0);
    //     BigInt temp;
    //     Xi[1] = bn_1;
    //     Xi[0] = -x[0];
    //     for(auto i = 2; i <= n; i++)
    //     {
    //         Xi[i] = Xi[i-1];
    //         for(auto j = i-1; j > 0; j--)
    //         {
    //             Xi[j] = (Xi[j] * (-x[i-1]) % order + Xi[j-1]) % order;
    //         }
    //         Xi[0] = Xi[0] * (-x[i-1]) % order;
    //     }
    //     for(auto i = 0; i < n; i++)
    //     {
    //         temp = bn_1;
    //         for(auto j = 0; j < n; j++)
    //         {
    //             if(i == j)
    //             {
    //                 continue;
    //             }
    //             temp = temp * (x[i] - x[j]) % order;
    //         }
    //         temp = temp.ModInverse(order);
    //         temp = y[i] * temp % order;

    //         for(auto j = 0; j <= n; j++)
    //         {
    //             Xi_temp[j] = Xi[j];
    //         }
    //         for(auto j = n-1; j >= 0; j--)
    //         {
    //             Fi[j] = Xi_temp[j+1];
    //             Xi_temp[j] = (Xi_temp[j] - (Xi_temp[j+1] * (-x[i])) % order) % order;
    //         }

    //         for(auto j = 0; j < n; j++)
    //         {
    //             result[j] = (result[j] + Fi[j] * temp %order) % order;
    //         }
            
    //     }
    //     return result;     
    // }
    /*
        Generate an argument PI for Relation 3 on pp.13: P = g^a h^b u^<a,b>
        transcript_str is introduced to be used as a sub-protocol
    */
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
            vec_C[d] = pp.vec_h[0] * P[0][d];
            for(auto j = 1; j < n; j++)
            {
                vec_C[d] += pp.vec_h[j] * P[j][d];
            }
            vec_C[d] += pp.h * rho[d];
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

        //compute P(x)
        std::vector<BigInt> vec_P_x(n);
        for(auto j = 0; j < n; j++)
        {
            vec_P_x[j] = BigInt(j);
        }
        std::vector<BigInt> vec_P = lagrange(vec_P_x, vec_s); 
        //std::vector<BigInt> vec_a(n);
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
        
        std::copy(vec_P.begin(), vec_P.end(), vec_y.begin());

        std::copy(vec_sk.begin(), vec_sk.end(), vec_y.begin()+n);
        std::copy(vec_r.begin(), vec_r.end(), vec_y.begin()+2*n);
        std::copy(vec_v.begin(), vec_v.end(), vec_y.begin()+3*n);

        BigInt rP = GenRandomBigIntLessThan(order);
        ECPoint P_equal = ECPointVectorMul(pp.vec_u, vec_y) + pp.u * rP;

        transcript_str = "";
        transcript_str = transcript_str + linbit_proof.A.ToByteString() + linbit_proof.C.ToByteString() + linbit_proof.D.ToByteString();
        // for(auto d = 0; d < m; d++)
        // {
        //     transcript_str += vec_C[d].ToByteString();
        // }
        // transcript_str += P_equal.ToByteString();

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
        zs = zs * x_pow.ModInverse(order) % order;
        zs = (-zs + order) % order;

        proof.P_equal = P_equal;
        
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

        AmorHom::PP amorhom_pp = AmorHom::Setup(n);
        amorhom_pp.g = pp.g;
        amorhom_pp.h = pp.h;
        amorhom_pp.u = pp.u;
        amorhom_pp.vec_h = pp.vec_h;
        amorhom_pp.vec_u = pp.vec_u;

        AmorHom::Instance amorhom_instance;
        amorhom_instance.P = P_equal;
        amorhom_instance.S = ECPointVectorMul(pp.vec_h, vec_p_eval) + right;
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

    /* Check if PI is a valid proof for eqmdl product statement (G1^w = H1 and G2^w = H2) */
    bool Verify(PP &pp, Instance &instance, std::string &transcript_str, Proof &proof)
    {
        bool Validity = false;
        transcript_str = "";
        size_t n = pp.VECTOR_LEN;
        size_t m = pp.LOG_VECTOR_LEN;
        //LinBit::PP linbit_pp = LinBit::Setup(pp.vec_g, pp.h, n);
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

        std::vector<BigInt> exp_x(m+1);
        exp_x[0] = bn_1;  
        for(auto k = 1; k <= m; k++){
            exp_x[k] = exp_x[k-1] * x % order; 
        }
        ECPoint right = proof.vec_C[0] * (bn_0 - exp_x[0] + order);
        for(auto d = 1;d < m; d++){
            right += proof.vec_C[d] * (bn_0 - exp_x[d] + order);
        }

        AmorHom::PP amorhom_pp = AmorHom::Setup(n);
        amorhom_pp.g = pp.g;
        amorhom_pp.h = pp.h;
        amorhom_pp.u = pp.u;
        amorhom_pp.vec_h = pp.vec_h;
        amorhom_pp.vec_u = pp.vec_u;

        AmorHom::Instance amorhom_instance;
        amorhom_instance.P = proof.P_equal;
        amorhom_instance.S = ECPointVectorMul(pp.vec_h, vec_p_eval) + right;
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
        bool V2 = AmorHom::Verify(amorhom_pp, amorhom_instance, transcript_str, amorhom_proof);

        Validity = V1 && V2;
        return Validity;
    }

}

#endif