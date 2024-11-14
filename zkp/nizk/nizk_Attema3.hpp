
#ifndef ATTEMA_3_HPP
#define ATTEMA_3_HPP

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"
#include "../nizk/nizk_lin_bit.hpp"
#include "../nizk/nizk_amorhom_for_Attema3.hpp"
#include "../nizk/nizk_kbit.hpp"
#include "../nizk/nizk_zkdl.hpp"

namespace Attema3
{

    using Serialization::operator<<;
    using Serialization::operator>>;

    size_t GetTheNthBit(size_t index, size_t n)
    {
        return (index>>n)&1;
    }

    std::vector<BigInt> poly_compute(std::vector<BigInt> &x, std::vector<BigInt> &y) 
    {
        size_t n = x.size();

        std::vector<std::vector<BigInt>> P;
        for(auto i=1;i<n;i++){
            std::vector<BigInt> A(2,bn_1);

            A[0] = -x[i];
            P.emplace_back(A);
        }

        std::vector<BigInt> res = PolyMul(P);
        
        BigInt temp = bn_1;
        for(auto i=1;i<n;i++){
            temp = temp * x[i].ModNegate(order) % order;
        }
        temp = temp.ModInverse(order);

        for(auto i=0;i<res.size();i++){
            res[i] = res[i] * temp % order;
        }

        return res;
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

    struct PP
    {
        size_t n;     // denotes the size of witness (witness is upto l = 2^VECTOR_LEN)
        size_t k;

        ECPoint g;
        ECPoint h;
        std::vector<ECPoint> vec_g;
    };

    struct Instance
    {   
        std::vector<ECPoint> vec_P;
    };

    struct Witness
    {
        std::vector<size_t> vec_S;
        std::vector<BigInt> vec_x;
    };

    struct Proof
    {
        ECPoint P;
        AmorHom::Proof amorhom_proof;
    };

    std::ofstream &operator<<(std::ofstream &fout, const Proof &proof)
    {
        fout << proof.P;
        fout << proof.amorhom_proof;
        return fout; 
    }

    std::ifstream &operator>>(std::ifstream &fin, Proof &proof)
    {
        fin >> proof.P;
        fin >> proof.amorhom_proof;
        return fin; 
    }

    PP Setup(size_t VECTOR_LEN, size_t K)
    {
        PP pp;

        pp.n = VECTOR_LEN;
        pp.k = K;
        pp.g = GenRandomGenerator();
        pp.h = GenRandomGenerator();
        pp.vec_g = GenRandomECPointVector(2*pp.n-pp.k);

        return pp;
    }


    Proof Prove(PP pp, Instance instance, Witness witness, std::string &transcript_str)
    {
        Proof proof;

        size_t n = pp.n;
        size_t k = pp.k;

        // prepare (x,y)
        std::vector<BigInt> vec_input(n-k+1, bn_0);
        std::vector<BigInt> vec_output(n-k+1, bn_0);

        vec_input[0] = bn_0;
        vec_output[0] = bn_1;
        size_t count = 1;
        for(auto i=1;i<=n;i++){
            bool flag = true;
            for(auto j=0;j<k;j++){
                if(witness.vec_S[j] == i){
                    flag = false;
                    break;
                }
            }
            if(flag == true){
                vec_input[count] = BigInt(i);
                count++;
            }
            if(count == n-k+1){
                break;
            }
        }

        // compute vec_a
        std::vector<BigInt> vec_a(n-k,bn_0);
        std::vector<BigInt> lagrange_res = poly_compute(vec_input, vec_output);

        for(auto i=0;i<n-k;i++){
            vec_a[i] = lagrange_res[i+1];
        }

        // compute P(x)
        std::vector<BigInt> vec_px(n+1,bn_0);
        vec_px[0] = bn_1;
        for(auto i=0;i<k;i++){
            size_t index = witness.vec_S[i];
            BigInt x_now = BigInt(index);
            std::vector<BigInt> vec_x_power = GenBigIntPowerVector(n-k+1,x_now);
            BigInt y_now = bn_1;
            for(auto j=1;j<=(n-k);j++){
                y_now = (y_now + lagrange_res[j] * vec_x_power[j] % order) % order;
            }
            vec_px[index] = y_now;
        }

        // set vec_y
        std::vector<BigInt> vec_y(2*n-k);
        std::copy(vec_a.begin(), vec_a.end(), vec_y.begin());
        for(auto i=0;i<n;i++){
            vec_y[i+n-k] = witness.vec_x[i] * vec_px[i+1]  % order;  
        }

        // pick tau
        BigInt tau = GenRandomBigIntLessThan(order);
        
        // compute P
        proof.P = ECPointVectorMul(pp.vec_g, vec_y) + pp.h * tau;

        // set init
        ECPoint init;
        init.SetInfinity();
        
        // run amorhom
        AmorHom::PP amorhom_pp = AmorHom::Setup(n);
        amorhom_pp.k = pp.k;
        amorhom_pp.g = pp.g;
        amorhom_pp.u = pp.h;
        amorhom_pp.vec_g.resize(2*n);  
        for(auto i=0;i<2*n-k;i++){
            amorhom_pp.vec_g[i] = pp.vec_g[i];
        }  
        for(auto i=2*n-k;i<2*n;i++){
            amorhom_pp.vec_g[i] = init;
        }

        AmorHom::Instance amorhom_instance;
        amorhom_instance.P = proof.P;
        amorhom_instance.vec_c = instance.vec_P;

        AmorHom::Witness amorhom_witness;
        amorhom_witness.rp = tau;
        amorhom_witness.vec_y.resize(2*n);
        for(auto i=0;i<2*n-k;i++){
            amorhom_witness.vec_y[i] = vec_y[i];
        }
        for(auto i=2*n-k;i<2*n;i++){
            amorhom_witness.vec_y[i] = bn_0;
        }

        transcript_str = "";
        AmorHom::Prove(amorhom_pp, amorhom_instance, amorhom_witness, transcript_str, proof.amorhom_proof);

        return proof;
    }

    bool Verify(PP &pp, Instance &instance, std::string &transcript_str, Proof &proof)
    {
        bool Validity = false;

        size_t n = pp.n;
        size_t k = pp.k;

        ECPoint init;
        init.SetInfinity();
        
        // run amorhom
        AmorHom::PP amorhom_pp = AmorHom::Setup(n);
        amorhom_pp.k = pp.k;
        amorhom_pp.g = pp.g;
        amorhom_pp.u = pp.h;

        amorhom_pp.vec_g.resize(2*n);
        for(auto i=0;i<2*n-k;i++){
            amorhom_pp.vec_g[i] = pp.vec_g[i];
        }
        for(auto i=2*n-k;i<2*n;i++){
            amorhom_pp.vec_g[i] = init;
        }

        AmorHom::Instance amorhom_instance;
        amorhom_instance.P = proof.P;
        amorhom_instance.vec_c = instance.vec_P;

        transcript_str = "";
        Validity = AmorHom::Verify(amorhom_pp, amorhom_instance, transcript_str, proof.amorhom_proof);

        return Validity;
    }

}

#endif
