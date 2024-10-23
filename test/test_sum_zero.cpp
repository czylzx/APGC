#define DEBUG

#include "../zkp/nizk/nizk_sum_zero.hpp"
#include "../crypto/setup.hpp"
#include "../commitment/pedersen.hpp"

void test_sum_zero(size_t N_max) {
    PrintSplitLine('-');  
    std::cout << "begin the test of sum_zero.hpp >>>" << std::endl; 

    // 初始化加密环境
    CRYPTO_Initialize(); 

    // 创建 Pedersen 承诺的公共参数
    Pedersen::PP pp_com = Pedersen::Setup(2);

    // 创建 SUM_ZERO 的公共参数
    SumZero::PP pp = SumZero::Setup(N_max, pp_com);

    // 创建实例和见证
    SumZero::Instance instance;
    SumZero::Witness witness;
    SumZero::Proof proof;



    // 生成见证值
    std::vector<BigInt> witness_values(N_max);
    

    // 生成承诺
    size_t t = 0;
    std::vector<BigInt> values(N_max, BigInt(t)); 

    instance.vec_com.push_back(Pedersen::Commit(pp_com, values, witness.values));
    

    // 生成证明
    std::string transcript_str = "";
    SumZero::Prove(pp, instance, witness, proof, transcript_str);

    // 验证证明
    transcript_str = "";
    bool testval = SumZero::Verify(pp, instance, proof, transcript_str);

    if (testval == true) {
        std::cout << "verify success" << std::endl;
    } else {
        std::cout << "verify fail" << std::endl;
    }

    // 清理加密环境
    CRYPTO_Finalize(); 
}

int main() {
    size_t N_max = 8;  // 最大数值数量
    test_sum_zero(N_max);

    return 0; 
}