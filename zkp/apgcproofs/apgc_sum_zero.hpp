#ifndef SUM_ZERO_HPP_
#define SUM_ZERO_HPP_

#include "../../crypto/ec_point.hpp"
#include "../../crypto/hash.hpp"
#include "../../commitment/pedersen.hpp"
#include "../../utility/polymul.hpp"
#include "../bulletproofs/innerproduct_proof.hpp"
#include <vector>
#include <iostream>

namespace SumZero
{

    using Serialization::operator<<;
    using Serialization::operator>>;

    struct PP
    {
        size_t n;                          
        ECPoint g;                      
    };

    struct Instance
    {
        ECPoint C; 
    };

    struct Witness
    {
        std::vector<BigInt> r;
        // BigInt r; 
    };

    struct Proof
    {
        ECPoint A;           
        BigInt z; 
    };

    PP Setup(size_t n)
    {
        PP pp;
        pp.n = n;
        pp.g = GenRandomGenerator();                                  
        return pp;
    }

    void Prove(PP &pp, Instance &instance, Witness &witness,  std::string &transcript_str,Proof &proof)
    {

        BigInt a = GenRandomBigIntLessThan(BigInt(order));
        proof.A = pp.g * a; 

        transcript_str += proof.A.ToByteString();

        BigInt x = Hash::StringToBigInt(transcript_str);

        BigInt sumr = bn_0;
        for(auto i=0;i<pp.n;i++)
        {
            sumr+= witness.r[i];
        }
        BigInt z = a.Add(x.Mul(sumr));
        // BigInt z = a.Add(x.Mul(witness.r));
        proof. z = z;

    }

    bool Verify(PP &pp, Instance &instance,  std::string &transcript_str,Proof &proof)
    {

        transcript_str = "";
        transcript_str += proof.A.ToByteString();
        BigInt x = Hash::StringToBigInt(transcript_str);
        ECPoint left = pp.g * proof.z;
        ECPoint right = proof.A + instance.C*x;
        if(left == right)
        {
            return true;
        } 
        else{
            return false;
        }
        


    } // namespace SumZero
}
#endif