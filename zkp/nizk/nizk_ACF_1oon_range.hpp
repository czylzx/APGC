
#ifndef ACF_1OON_RANGE_HPP
#define ACF_1OON_RANGE_HPP

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"
#include "../nizk/nizk_amorhom_for_ACF1oonRange.hpp"
#include "../nizk/nizk_kbit.hpp"
#include "../nizk/nizk_zkdl.hpp"
#include "../bulletproofs/bullet_proof.hpp"

namespace ACF_1oon_range
{

    using Serialization::operator<<;
    using Serialization::operator>>;

    size_t GetTheNthBit(size_t index, size_t n)
    {
        return (index>>n)&1;
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
        size_t VECTOR_LEN;   
        ECPoint g;
        ECPoint h;
        ECPoint u;
        ECPoint u_new;
        std::vector<ECPoint> vec_g;
        std::vector<ECPoint> vec_h;
        std::vector<ECPoint> vec_u;
        // for innerproduct in bullet
        std::vector<ECPoint> vec_ip_g;
        std::vector<ECPoint> vec_ip_h;
    };

    struct Instance
    {   
        std::vector<ECPoint> vec_P;
    };

    struct Witness
    {
        size_t l;
        std::vector<BigInt> vec_x;
    };

    struct Proof
    {
        ECPoint P;
        ECPoint S;
        ECPoint V;
        std::vector<ECPoint> vec_C;
        BigInt zs;
        AmorHom::Proof amorhom_proof;
        Kbit::Proof kbit_proof;
        BigInt k;
        ZkdlProduct::Proof zkdl_proof;
        BigInt f_zkdl;
        ECPoint Ax;
        Bullet::Proof bullet_proof;
    }; 

    std::ofstream &operator<<(std::ofstream &fout, const Proof &proof)
    {
        fout << proof.P << proof.S << proof.vec_C << proof.zs << proof.k << proof.f_zkdl << proof.Ax; 
        fout << proof.amorhom_proof;   
        fout << proof.kbit_proof;   
        fout << proof.zkdl_proof;   
        return fout;
    }

    std::ifstream &operator>>(std::ifstream &fin, Proof &proof)
    {
        fin >> proof.P >> proof.S >> proof.vec_C >> proof.zs >> proof.k >> proof.f_zkdl >> proof.Ax; 
        fin >> proof.amorhom_proof;   
        fin >> proof.kbit_proof;   
        fin >> proof.zkdl_proof;   
        return fin;
    }     

    PP Setup(size_t VECTOR_LEN)
    {
        PP pp;
        pp.VECTOR_LEN = VECTOR_LEN;
        pp.g = GenRandomGenerator();
        pp.h = GenRandomGenerator();
        pp.u = GenRandomGenerator();
        pp.u_new = GenRandomGenerator();
        pp.vec_g = GenRandomECPointVector(VECTOR_LEN);
        pp.vec_h = GenRandomECPointVector(VECTOR_LEN);
        pp.vec_u = GenRandomECPointVector(VECTOR_LEN*3);
        pp.vec_ip_g = GenRandomECPointVector(32);
        pp.vec_ip_h = GenRandomECPointVector(32);
        return pp;
    }


    Proof Prove(PP pp, Instance instance, Witness witness, std::string &transcript_str)
    {
        Proof proof;

        size_t n = pp.VECTOR_LEN;

        // set vec_s
        std::vector<BigInt> vec_s(n,bn_0);
        vec_s[witness.l] = bn_1;

        // pick rP,rS,rV
        BigInt rP = GenRandomBigIntLessThan(order);
        BigInt rS = GenRandomBigIntLessThan(order);
        BigInt rV = GenRandomBigIntLessThan(order);

        // prepare for set vec_y
        std::vector<BigInt> vec_sx(2*n,bn_0);
        vec_sx[2*witness.l] = witness.vec_x[0];
        vec_sx[2*witness.l+1] = witness.vec_x[1];

        std::vector<BigInt> vec_srV(n,bn_0);
        vec_srV[witness.l] = rV;

        // set vec_y
        std::vector<BigInt> vec_y(4*n,bn_0);
        std::copy(vec_s.begin(), vec_s.end(), vec_y.begin());
        std::copy(vec_sx.begin(), vec_sx.end(),vec_y.begin()+n);
        std::copy(vec_srV.begin(), vec_srV.end(),vec_y.begin()+3*n);

        // prepare for compute
        std::vector<ECPoint> A;
        std::vector<BigInt> B;

        // compute V
        proof.V = pp.g * rV + pp.h * witness.vec_x[1];      

        // compute P
        A.resize(4*n+1);
        B.resize(4*n+1);

        std::copy(pp.vec_g.begin(), pp.vec_g.end(), A.begin());
        std::copy(pp.vec_u.begin(), pp.vec_u.end(), A.begin()+n);
        A[4*n] = pp.u;

        std::copy(vec_y.begin(), vec_y.end(), B.begin());
        B[4*n] = rP;

        proof.P = ECPointVectorMul(A, B);

        // compute vec_s - vec_1
        std::vector<BigInt> vec_s_1(n, (-bn_1));
        vec_s_1[witness.l] = bn_0;

        // compute S
        A.resize(2*n+1);
        B.resize(2*n+1);

        std::copy(pp.vec_g.begin(), pp.vec_g.end(), A.begin());
        std::copy(pp.vec_h.begin(), pp.vec_h.end(), A.begin()+n);
        A[2*n] = pp.u;
        
        std::copy(vec_s.begin(), vec_s.end(), B.begin());
        std::copy(vec_s_1.begin(), vec_s_1.end(), B.begin()+n);
        B[2*n] = rS;  

        proof.S = ECPointVectorMul(A, B);

        // range proof
        size_t range_len = 32;
        size_t max_agg_num = 1;
        
        // Bullet proof v
        Bullet::PP Bullet_pp = Bullet::Setup(range_len,max_agg_num);
        Bullet::Instance Bullet_instance;
        Bullet::Witness Bullet_witness;

        // set pp
        Bullet_pp.g = pp.g;
        Bullet_pp.h = pp.h;
        Bullet_pp.u = pp.u;
        Bullet_pp.vec_g = pp.vec_ip_g;
        Bullet_pp.vec_h = pp.vec_ip_h;
        
        // set instance
        Bullet_instance.C.resize(1);
        Bullet_instance.C[0] = proof.V;

        // set witness
        Bullet_witness.r = GenRandomBigIntVectorLessThan(max_agg_num,order);
        Bullet_witness.v = GenRandomBigIntVectorLessThan(max_agg_num,order);
        Bullet_witness.r[0] = rV;
        Bullet_witness.v[0] = witness.vec_x[1];

        // set transcript_str
        std::string Bullet_str = "";
        
        // call Bullet proof
        Bullet::Prove(Bullet_pp, Bullet_instance, Bullet_witness, Bullet_str, proof.bullet_proof);

        // run kbit
        Kbit::PP kbit_pp = Kbit::Setup(n);
        kbit_pp.g = pp.g;
        kbit_pp.h = pp.h;
        kbit_pp.u = pp.u;
        kbit_pp.u_new = pp.u_new;
        kbit_pp.vec_h = pp.vec_h;
        kbit_pp.vec_g = pp.vec_g;

        Kbit::Instance kbit_instance;
        kbit_instance.P = proof.S;
        kbit_instance.k = 1;

        Kbit::Witness kbit_witness;
        kbit_witness.vec_a = vec_s;
        kbit_witness.vec_b = vec_s_1;
        kbit_witness.r = rS;

        transcript_str = "";
        proof.kbit_proof = Kbit::Prove(kbit_pp, kbit_instance, kbit_witness, transcript_str);

        // run zkdls
        // set pp
        ZkdlProduct::PP zkdl_pp = ZkdlProduct::Setup(4*n, false);
        zkdl_pp.vec_g.resize(4*n);
        std::copy(pp.vec_h.begin(), pp.vec_h.end(), zkdl_pp.vec_g.begin());
        std::copy(pp.vec_u.begin(), pp.vec_u.end(), zkdl_pp.vec_g.begin()+n);

        // pick random
        std::vector<BigInt> a_init = GenRandomBigIntVectorLessThan(4*n, order);
        BigInt b_init = GenRandomBigIntLessThan(order);

        // compute A
        proof.Ax = ECPointVectorMul(zkdl_pp.vec_g, a_init) + pp.u * b_init;

        // compute e
        BigInt e_zkdl = Hash::StringToBigInt(proof.Ax.ToByteString());

        // compute f
        proof.f_zkdl = (b_init + e_zkdl * (rP - rS)) % order;

        // prepare vec_a
        std::vector<BigInt> vec_a_zkdl(4*n);
        for(auto i=0;i<n;i++){
            vec_a_zkdl[i] = (bn_1 - vec_s[i]) % order;
        }
        std::copy(vec_y.begin()+n, vec_y.end(), vec_a_zkdl.begin()+n);
        std::vector<BigInt> vec_y_zkdl_temp = BigIntVectorScalar(vec_a_zkdl, e_zkdl);
        std::vector<BigInt> vec_y_zkdl = BigIntVectorModAdd(a_init, vec_y_zkdl_temp, order);

        ZkdlProduct::Instance zkdl_instance;
        zkdl_instance.G = (proof.P + proof.S.Invert()) * e_zkdl + proof.Ax - pp.u * proof.f_zkdl;
        
        ZkdlProduct::Witness zkdl_witness;
        zkdl_witness.vec_a = vec_y_zkdl;

        transcript_str = "";
        ZkdlProduct::Prove(zkdl_pp, zkdl_instance,zkdl_witness, transcript_str, proof.zkdl_proof);

        // run amorhom
        std::vector<ECPoint> vec_g_amorhom(4*n);
        std::copy(pp.vec_g.begin(), pp.vec_g.end(), vec_g_amorhom.begin());
        std::copy(pp.vec_u.begin(), pp.vec_u.end(), vec_g_amorhom.begin()+n);

        AmorHom::PP amorhom_pp = AmorHom::Setup(n);
        amorhom_pp.g = pp.g;
        amorhom_pp.h = pp.h;
        amorhom_pp.u = pp.u;
        amorhom_pp.vec_g = vec_g_amorhom;

        AmorHom::Instance amorhom_instance;
        amorhom_instance.P = proof.P;
        amorhom_instance.V = proof.V;
        amorhom_instance.vec_c = instance.vec_P;

        AmorHom::Witness amorhom_witness;
        amorhom_witness.vec_y = vec_y;
        amorhom_witness.rp = rP;
        transcript_str = "";
        AmorHom::Prove(amorhom_pp, amorhom_instance, amorhom_witness, transcript_str, proof.amorhom_proof);

        #ifdef DEBUG
        std::cerr << "ACF Koon Proof Generation Finishes >>>" << std::endl;
        #endif 

        return proof;
    }

    bool Verify(PP &pp, Instance &instance, std::string &transcript_str, Proof &proof)
    {
        bool Validity = false;
        std::vector<bool> vec_condition(4, false);

        size_t n = pp.VECTOR_LEN;

        // run range
        // range proof
        size_t range_len = 32;
        size_t max_agg_num = 1;
        
        // Bullet proof v
        Bullet::PP Bullet_pp = Bullet::Setup(range_len,max_agg_num);
        Bullet::Instance Bullet_instance;
        Bullet::Witness Bullet_witness;

        // set pp
        Bullet_pp.g = pp.g;
        Bullet_pp.h = pp.h;
        Bullet_pp.u = pp.u;
        Bullet_pp.vec_g = pp.vec_ip_g;
        Bullet_pp.vec_h = pp.vec_ip_h;
        
        // set instance
        Bullet_instance.C.resize(1);
        Bullet_instance.C[0] = proof.V;

        // set transcript_str
        std::string Bullet_str = "";
        
        // call Bullet proof
        vec_condition[0] = Bullet::Verify(Bullet_pp, Bullet_instance, Bullet_str, proof.bullet_proof);
    
        // run kbit
        Kbit::PP kbit_pp = Kbit::Setup(n);
        kbit_pp.g = pp.g;
        kbit_pp.h = pp.h;
        kbit_pp.u = pp.u;
        kbit_pp.u_new = pp.u_new;
        kbit_pp.vec_h = pp.vec_h;
        kbit_pp.vec_g = pp.vec_g;

        Kbit::Instance kbit_instance;
        kbit_instance.P = proof.S;
        kbit_instance.k = 1;

        transcript_str = "";
        vec_condition[1] = Kbit::Verify(kbit_pp, kbit_instance, transcript_str, proof.kbit_proof);

        // run zkdl
        ZkdlProduct::PP zkdl_pp = ZkdlProduct::Setup(4*n, false);
        zkdl_pp.vec_g.resize(4*n);
        std::copy(pp.vec_h.begin(), pp.vec_h.end(), zkdl_pp.vec_g.begin());
        std::copy(pp.vec_u.begin(), pp.vec_u.end(), zkdl_pp.vec_g.begin()+n);

        ZkdlProduct::Instance zkdl_instance;
        BigInt e_zkdl = Hash::StringToBigInt(proof.Ax.ToByteString());

        zkdl_instance.G = (proof.P + proof.S.Invert()) * e_zkdl + proof.Ax - pp.u * proof.f_zkdl;

        transcript_str = "";
        vec_condition[2] = ZkdlProduct::Verify(zkdl_pp, zkdl_instance, transcript_str, proof.zkdl_proof);

        // run amorhom
        std::vector<ECPoint> vec_g_amorhom(4*n);
        std::copy(pp.vec_g.begin(), pp.vec_g.end(), vec_g_amorhom.begin());
        std::copy(pp.vec_u.begin(), pp.vec_u.end(), vec_g_amorhom.begin()+n);

        AmorHom::PP amorhom_pp = AmorHom::Setup(n);
        amorhom_pp.g = pp.g;
        amorhom_pp.h = pp.h;
        amorhom_pp.u = pp.u;
        amorhom_pp.vec_g = vec_g_amorhom;


        AmorHom::Instance amorhom_instance;
        amorhom_instance.P = proof.P;
        amorhom_instance.V = proof.V;
        amorhom_instance.vec_c = instance.vec_P;

        transcript_str = "";
        vec_condition[3] = AmorHom::Verify(amorhom_pp, amorhom_instance, transcript_str, proof.amorhom_proof);

        Validity = vec_condition[0] && vec_condition[1] && vec_condition[2] && vec_condition[3];



        #ifdef DEBUG
        for(auto i = 0; i < 4; i++){
            std::cout << std::boolalpha << "Condition "<< std::to_string(i) <<" (ACF 1oon range proof) = " 
                  << vec_condition[i] << std::endl; 
        }

        if (Validity){ 
            std::cout << "NIZK proof for ACF 1oon range accepts >>>" << std::endl; 
        } else {
            std::cout << "NIZK proof for ACF 1oon range rejects >>>" << std::endl; 
        }
        #endif

        return Validity;
    }

}

#endif
