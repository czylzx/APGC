#ifndef AmorHom_HPP_
#define AmorHom_HPP_

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"
#include "../../commitment/pedersen.hpp"
#include "../../utility/polymul.hpp"
#include "../bulletproofs/innerproduct_proof.hpp"
#include "./apgc_eqmdl4two_proof.hpp"
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
        std::vector<ECPoint> vec_h;
        std::vector<ECPoint> vec_u;                    
    };

    struct Instance
    {
        ECPoint P; 
        ECPoint S;
        std::vector<ECPoint> vec_pk;
        std::vector<ECPoint> Sum_CL;
        std::vector<ECPoint> Sum_CR;
        ECPoint Refresh_CL;
        ECPoint Refresh_CR;
        BigInt x;
        BigInt zs;
    };

    struct Witness
    {
        std::vector<BigInt> vec_y;
        BigInt rp;
        // BigInt r; 
    };

    struct Proof
    {
        EqmdlProduct2::Proof eqmdl_proof;
        ECPoint Ap;
        ECPoint Af;
        BigInt f;
    };

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
        transcript_str += instance.Refresh_CL.ToByteString();
        transcript_str += instance.Refresh_CR.ToByteString();
        for(auto i = 0; i < pp.n; i++){
            transcript_str += instance.vec_pk[i].ToByteString();
        }

        for(auto i = 0; i < pp.n; i++){
            transcript_str += instance.Sum_CL[i].ToByteString();
        }

        BigInt p = Hash::StringToBigInt(transcript_str);

        //std::vector<ECPoint> vec_base_new(pp.n*4);

        ECPoint init;
        init.SetInfinity();
        
        std::vector<ECPoint> vec_base_new(pp.n*4, init);

        std::vector<std::vector<ECPoint>> vec_y_base_1(pp.n, std::vector<ECPoint>(pp.n*4,init));
        std::vector<std::vector<ECPoint>> vec_y_base_2(pp.n, std::vector<ECPoint>(pp.n*4,init));
        std::vector<std::vector<ECPoint>> vec_y_base_3(pp.n, std::vector<ECPoint>(pp.n*4,init));
        std::vector<std::vector<ECPoint>> vec_y_base_4(pp.n, std::vector<ECPoint>(pp.n*4,init));
        std::vector<ECPoint> vec_y_base_5(pp.n*4,init);
        

        size_t n = pp.n;
        

        //std::vector<BigInt> vec_j_power = GenBigIntPowerVector(pp.n, j);
        std::vector<std::vector<BigInt>> vec_j_power(pp.n, std::vector<BigInt>(pp.n));
        for(auto j = 0; j < pp.n; j++)
        {
            if(j == 0)
            {
                for(auto i = 0; i < pp.n; i++)
                {
                    vec_j_power[j][i] = bn_1;
                }
            }
            else
            {
                vec_j_power[j] = GenBigIntPowerVector(pp.n, j);
            }
            
        }
        std::vector<BigInt> vec_p = GenBigIntPowerVector(pp.n*4+1, p);
        auto start_time = std::chrono::steady_clock::now();
        for(auto j=0; j < pp.n; j++)
        {
            vec_y_base_1[j][j+n] =  vec_y_base_1[j][j+n] + pp.g;
            
            vec_y_base_2[j][j+2*n] =  vec_y_base_2[j][j+2*n] + instance.vec_pk[j];
            
            vec_y_base_3[j][j+2*n] =  vec_y_base_3[j][j+2*n] + pp.g;
            vec_y_base_3[j][j+3*n] =  vec_y_base_3[j][j+3*n] + pp.h;
            

            vec_y_base_4[j][j+n] =  vec_y_base_4[j][j+n] + instance.Refresh_CR - instance.Sum_CR[j];
           

            for(auto i = 0; i < pp.n; i++)
            {
               vec_y_base_1[j][i] =  vec_y_base_1[j][i] + instance.vec_pk[j] * (-vec_j_power[j][i]);
               vec_y_base_2[j][i] =  vec_y_base_2[j][i] + instance.Refresh_CL * (-vec_j_power[j][i]);
               vec_y_base_3[j][i] =  vec_y_base_3[j][i] + instance.Refresh_CR * (-vec_j_power[j][i]);
               vec_y_base_4[j][i] =  vec_y_base_4[j][i] + (instance.Refresh_CL - instance.Sum_CL[j]) * (-vec_j_power[j][i]);
               vec_y_base_5[i] =  vec_y_base_5[i] + pp.vec_h[j] * (vec_j_power[j][i]);
            }
           
        }
        auto end_time = std::chrono::steady_clock::now();
        auto time_diff = end_time - start_time;
        std::cout << "time for computing vec_y_base = " << std::chrono::duration<double, std::milli>(time_diff).count() << " ms" << std::endl;
        // std::vector<ECPoint> vec_base_new_text(pp.n*4, init);
        start_time = std::chrono::steady_clock::now();
        for(auto i = 0; i < pp.n*4; i++)
        {
            for(auto j = 0; j < pp.n; j++)
            {
               vec_base_new[i] = vec_base_new[i] + vec_y_base_1[j][i] * vec_p[j]
                                + vec_y_base_2[j][i] * vec_p[j+n]
                                + vec_y_base_3[j][i] * vec_p[j+2*n]
                                + vec_y_base_4[j][i] * vec_p[j+3*n];
            //   vec_base_new_text[i] = vec_base_new_text[i] + vec_y_base_1[j][i] * vec_p[j]
            //                     + vec_y_base_2[j][i] * vec_p[j+n]
            //                     + vec_y_base_3[j][i] * vec_p[j+2*n]
            //                     + vec_y_base_4[j][i] * vec_p[j+3*n];              
            }
            vec_base_new[i] = vec_base_new[i] + vec_y_base_5[i] * vec_p[pp.n*4];    
        }
        end_time = std::chrono::steady_clock::now();
        time_diff = end_time - start_time;
        std::cout << "time for computing vec_base_new = " << std::chrono::duration<double, std::milli>(time_diff).count() << " ms" << std::endl;
        // ECPoint test_F = ECPointVectorMul(vec_base_new_text, witness.vec_y); 
        // PrintECPointVector(vec_base_new_text, "vec_base_new_text");
        // PrintBigIntVector(witness.vec_y, "witness.vec_y");
        // if(test_F.IsAtInfinity())
        // {
        //     std::cout << "test_F is at infinity" << std::endl;
        // }
        // else
        // {
        //     std::cout << "test_F is not at infinity" << std::endl;
        // }
        //compute F'
        ECPoint F_prime;
        BigInt x_exp_m = instance.x.ModExp(-pp.m, order);
        F_prime = instance.S * x_exp_m + pp.h * (instance.zs);
        F_prime = F_prime * vec_p[pp.n*4];

        // ECPoint F_test = ECPointVectorMul(vec_base_new, witness.vec_y);
        
        // if(F_test != F_prime)
        // {
        //     std::cout << "F_test != F_prime" << std::endl;
        // }
        // else
        // {
        //     std::cout << "F_test == F_prime" << std::endl;
        // }

        F_prime = ECPointVectorMul(vec_base_new, witness.vec_y);
        std::vector<BigInt> vec_a_init = GenRandomBigIntVectorLessThan(pp.n*4, order);
        BigInt b_inint = GenRandomBigIntLessThan(order);

        ECPoint Ap = ECPointVectorMul(pp.vec_u, vec_a_init) +pp.u * b_inint;
        ECPoint Af = ECPointVectorMul(vec_base_new, vec_a_init);

        proof.Af = Af;
        proof.Ap = Ap;

        transcript_str = "";
        transcript_str += Ap.ToByteString() + Af.ToByteString();

        //compute the challenge e
        BigInt e = Hash::StringToBigInt(transcript_str);
        
        std::vector<BigInt> vec_z_tmp =BigIntVectorModScalar(witness.vec_y, e, order);

        std::vector<BigInt> vec_z = BigIntVectorModAdd(vec_z_tmp, vec_a_init, order);

        BigInt f = (b_inint + e * witness.rp) % order;

        proof.f = f;
        //std::vector<BigInt> vec_z = BigIntVectorModScalar(witness.vec_y, e, order);     
        EqmdlProduct2:: PP eqmdl_pp = EqmdlProduct2::Setup(pp.n*4, true);
        eqmdl_pp.vec_g = vec_base_new;
        eqmdl_pp.vec_p = pp.vec_u;

        //PrintECPointVector(eqmdl_pp.vec_g, "eqmdl_pp.vec_g");
        //PrintECPointVector(eqmdl_pp.vec_p, "eqmdl_pp.vec_p");
        EqmdlProduct2::Instance eqmdl_instance;
        eqmdl_instance.P = instance.P * e + Ap + pp.u * (-f);
        eqmdl_instance.G = F_prime * e + Af;

        /*text*/
        // ECPoint left1 = ECPointVectorMul(pp.vec_u, vec_z);
        // ECPoint left2 = ECPointVectorMul(vec_base_new, vec_z);

        // if(left1 != eqmdl_instance.P)
        // {
        //     std::cout << "left1 != eqmdl_instance.P" << std::endl;
        // }
        // if(left2 != eqmdl_instance.G)
        // {
        //     std::cout << "left2 != eqmdl_instance.G" << std::endl;
        // }
        //eqmdl_instance.P.Print("eqmdl_instance.P");
        //eqmdl_instance.G.Print("eqmdl_instance.G");

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
        transcript_str += instance.Refresh_CL.ToByteString();
        transcript_str += instance.Refresh_CR.ToByteString();
        for(auto i = 0; i < pp.n; i++){
            transcript_str += instance.vec_pk[i].ToByteString();
        }
     
        for(auto i = 0; i < pp.n; i++){
            transcript_str += instance.Sum_CL[i].ToByteString();
        }

        BigInt p = Hash::StringToBigInt(transcript_str);

        //std::vector<ECPoint> vec_base_new(pp.n*4);

        ECPoint init;
        init.SetInfinity();
        std::vector<ECPoint> vec_base_new(pp.n*4, init);

        std::vector<std::vector<ECPoint>> vec_y_base_1(pp.n, std::vector<ECPoint>(pp.n*4,init));
        std::vector<std::vector<ECPoint>> vec_y_base_2(pp.n, std::vector<ECPoint>(pp.n*4,init));
        std::vector<std::vector<ECPoint>> vec_y_base_3(pp.n, std::vector<ECPoint>(pp.n*4,init));
        std::vector<std::vector<ECPoint>> vec_y_base_4(pp.n, std::vector<ECPoint>(pp.n*4,init));
        std::vector<ECPoint> vec_y_base_5(pp.n*4,init);
        

        size_t n = pp.n;

        //std::vector<BigInt> vec_j_power = GenBigIntPowerVector(pp.n, j);
        std::vector<std::vector<BigInt>> vec_j_power(pp.n, std::vector<BigInt>(pp.n));
        for(auto j = 0; j < pp.n; j++)
        {
            if(j == 0)
            {
                for(auto i = 0; i < pp.n; i++)
                {
                    vec_j_power[j][i] = bn_1;
                }
            }
            else
            {
                vec_j_power[j] = GenBigIntPowerVector(pp.n, j);
            }
            
        }
        std::vector<BigInt> vec_p = GenBigIntPowerVector(pp.n*4+1, p);
  
        for(auto j=0; j < pp.n; j++)
        {
            vec_y_base_1[j][j+n] =  vec_y_base_1[j][j+n] + pp.g;
            
            vec_y_base_2[j][j+2*n] =  vec_y_base_2[j][j+2*n] + instance.vec_pk[j];
            
            vec_y_base_3[j][j+2*n] =  vec_y_base_3[j][j+2*n] + pp.g;
            vec_y_base_3[j][j+3*n] =  vec_y_base_3[j][j+3*n] + pp.h;
            

            vec_y_base_4[j][j+n] =  vec_y_base_4[j][j+n] + instance.Refresh_CR - instance.Sum_CR[j];
           

            for(auto i = 0; i < pp.n; i++)
            {
               vec_y_base_1[j][i] =  vec_y_base_1[j][i] + instance.vec_pk[j] * (-vec_j_power[j][i]);
               vec_y_base_2[j][i] =  vec_y_base_2[j][i] + instance.Refresh_CL * (-vec_j_power[j][i]);
               vec_y_base_3[j][i] =  vec_y_base_3[j][i] + instance.Refresh_CR * (-vec_j_power[j][i]);
               vec_y_base_4[j][i] =  vec_y_base_4[j][i] + (instance.Refresh_CL - instance.Sum_CL[j]) * (-vec_j_power[j][i]);
               vec_y_base_5[i] =  vec_y_base_5[i] + pp.vec_h[j] * (vec_j_power[j][i]);
            }
           
        }
        for(auto i = 0; i < pp.n*4; i++)
        {
            for(auto j = 0; j < pp.n; j++)
            {
               vec_base_new[i] = vec_base_new[i] + vec_y_base_1[j][i] * vec_p[j]
                                + vec_y_base_2[j][i] * vec_p[j+n]
                                + vec_y_base_3[j][i] * vec_p[j+2*n]
                                + vec_y_base_4[j][i] * vec_p[j+3*n];
                                
            }
            vec_base_new[i] = vec_base_new[i] + vec_y_base_5[i] * vec_p[pp.n*4];    
        }
    
        //compute F'
        ECPoint F_prime;
        BigInt x_exp_m = instance.x.ModExp(-pp.m, order);
        F_prime = instance.S * x_exp_m + pp.h * (instance.zs);
        F_prime = F_prime * vec_p[pp.n*4];

        transcript_str = "";
        transcript_str += proof.Ap.ToByteString() + proof.Af.ToByteString();

        //compute the challenge e
        BigInt e = Hash::StringToBigInt(transcript_str);

        EqmdlProduct2:: PP eqmdl_pp = EqmdlProduct2::Setup(pp.n*4, true);
        eqmdl_pp.vec_g = vec_base_new;
        eqmdl_pp.vec_p = pp.vec_u;

        //PrintECPointVector(eqmdl_pp.vec_g, "vec_g");
        //PrintECPointVector(eqmdl_pp.vec_p, "vec_p");
   
        EqmdlProduct2::Instance eqmdl_instance;
        eqmdl_instance.P = instance.P * e + proof.Ap + pp.u * (-proof.f);
        eqmdl_instance.G = F_prime * e + proof.Af ;
    
        //eqmdl_instance.P.Print("P");
        //eqmdl_instance.G.Print("G");
        EqmdlProduct2::Proof eqmdl_proof;

        transcript_str = "";
        eqmdl_proof = proof.eqmdl_proof;
        bool Validity = EqmdlProduct2::Verify(eqmdl_pp, eqmdl_instance, transcript_str, eqmdl_proof);

        // if(Validity == false)
        // {
        //    std::cout << "amorhom_proof is invalid" << std::endl;
        // }
        // else
        // {
        //     std::cout << "amorhom_proof is valid" << std::endl;
        // }
        return Validity;

    } 
}
#endif