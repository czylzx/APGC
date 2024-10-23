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

    // 定义公共参数结构，包含椭圆曲线生成元和其他必要参数
    struct PP
    {
        size_t n;                          // 数组长度，表示见证值的数量
        Pedersen::PP com_part;             // Pedersen承诺的公共参数
        ECPoint g, h;                      // 椭圆曲线生成元
        std::vector<ECPoint> vec_g, vec_h; // 用于内积证明的公共参数
    };

    // 实例结构，包含一系列承诺
    struct Instance
    {
        std::vector<ECPoint> vec_com; // 承诺向量
    };

    // 见证结构，包含需要证明的数值
    struct Witness
    {
        BigInt values; // 需要证明的数值
    };

    // 证明结构，包含证明过程中生成的承诺和响应
    struct Proof
    {
        ECPoint A;           // 承诺A
        BigInt taux, mu, tx; // 挑战响应
        // InnerProduct::Proof ip_proof; // 内积证明
    };

    // Setup阶段：初始化公共参数
    PP Setup(size_t n, Pedersen::PP &com_part)
    {
        PP pp;
        pp.n = n;
        pp.com_part = com_part;
        pp.g = generator;                                  // 选择一个生成元g
        pp.h = Hash::StringToECPoint(pp.g.ToByteString()); // 通过哈希生成另一个生成元h
        pp.vec_g.resize(n, pp.g);                          // 初始化vec_g
        pp.vec_h.resize(n, pp.h);                          // 初始化vec_h
        return pp;
    }

    // Prove阶段：证明者生成证明
    void Prove(PP &pp, Instance &instance, Witness &witness, Proof &proof, std::string &transcript_str)
    {

        // 随机选择a，用于生成承诺A
        BigInt a = GenRandomBigIntLessThan(BigInt(order));
        proof.A = pp.g * a; // 计算承诺A = g^a

        // 更新transcript字符串，加入承诺A
        transcript_str += proof.A.ToByteString();

        // 计算挑战e
        BigInt x = Hash::StringToBigInt(transcript_str);

        // proof.z = (a + e * witness.values) % order; // z = a+e*w mod q
        // 生成响应z
        BigInt r = witness.values;
        BigInt z = (a + x.Mul(r)) % order;
        proof.taux = z; // 设置响应taux
    }

    // Verify阶段：验证者验证证明
    bool Verify(PP &pp, Instance &instance, Proof &proof, std::string &transcript_str)
    {

        // 重新计算挑战e
        transcript_str += proof.A.ToByteString();
        BigInt e = Hash::StringToBigInt(transcript_str);

        // 计算预期的承诺A
        ECPoint A_expected = pp.g * proof.taux;

        // 初始化 InnerProduct::PP 和 InnerProduct::Instance 对象
        // InnerProduct::PP ip_pp = InnerProduct::Setup(pp.n, false);
        // InnerProduct::Instance ip_instance;

        bool validity = proof.A == A_expected;

        return validity;

    } // namespace SumZero
}
#endif