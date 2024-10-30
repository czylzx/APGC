/****************************************************************************
this hpp implements the APGC functionality 
*****************************************************************************/
#ifndef APGC_HPP_
#define APGC_HPP_

#include "../pke/twisted_exponential_elgamal.hpp"        // implement Twisted ElGamal  
#include "../zkp/nizk/nizk_plaintext_equality.hpp" // NIZKPoK for plaintext equality
#include "../zkp/nizk/nizk_plaintext_knowledge.hpp"        // NIZKPoK for ciphertext/honest encryption 
#include "../zkp/nizk/nizk_dlog_equality.hpp"      // NIZKPoK for dlog equality
#include "../zkp/nizk/nizk_dlog_knowledge.hpp"     // NIZKPoK for dlog knowledge
#include "../zkp/bulletproofs/bullet_proof.hpp"    // implement Log Size Bulletproof
#include "../gadget/range_proof.hpp"
#include "../utility/serialization.hpp"
#include "../zkp/apgcproofs/apgc_wellform.hpp"
#include "../zkp/nizk/nizk_multi_plaintext_equality.hpp"
#include "../zkp/apgcproofs/apgc_sum_zero.hpp"
#include "../zkp/nizk/nizk_any_out_of_many.hpp"
#include "../zkp/apgcproofs/apgc_sdr_solvent_equal.hpp"

#define DEMO           // demo mode 
//#define DEBUG        // show debug information 


//!!!!!!!!!!!!!!!! if you want to invoke the Bulletproofs, must note taux and tx's order ! The order is reversed !!

namespace APGC{

using Serialization::operator<<; 
using Serialization::operator>>; 

// define the structure of system parameters

struct PP{    
    size_t SN_LEN;    // sn length
    size_t MAX_RECEIVER_NUM; // number of maximum receivers (for now, we require this value to be 2^n - 1)
    BigInt MAXIMUM_COINS; 
    size_t SLACK_PARTICIPANT_NUM; // number of slack participants

    Bullet::PP bullet_part; 
    TwistedExponentialElGamal::PP enc_part;

    ECPoint pka; // supervisor's pk
};

// define the structure of system parameters
struct SP{
    BigInt ska;   // supervisor's sk
};

struct Account{
    std::string identity;     // id
    ECPoint pk;              // public key
    BigInt sk;              // secret key
    TwistedExponentialElGamal::CT balance_ct;  // current balance
    BigInt m;               // dangerous (should only be used for speeding up the proof generation)
    //BigInt sn; 
};

struct SuperviseResult{
    ECPoint sender_pk;
    std::vector<ECPoint> receiver_pks;
    std::vector<BigInt> receiver_coins_values;
};

template <typename T>
std::string GetCTxFileName(T &newCTx)
{
    std::string ctx_file = newCTx.sn.ToHexString()+".ctx"; 
    return ctx_file; 
}

void PrintPP(PP &pp)
{
    PrintSplitLine('-');
    std::cout << "pp content >>>>>>" << std::endl; 
    std::cout << "MAX_RECEIVER_NUM = " << pp.MAX_RECEIVER_NUM << std::endl; // number of sub-argument (for now, we require m to be the power of 2)
    std::cout << "SN_LEN = " << pp.SN_LEN << std::endl; 

    pp.pka.Print("supervisor's pk"); 
    
    PrintSplitLine('-'); 
}

void PrintAccount(Account &Acct)
{
    std::cout << Acct.identity << " account information >>> " << std::endl;     
    Acct.pk.Print("pk"); 
    std::cout << "encrypted balance:" << std::endl; 
    TwistedExponentialElGamal::PrintCT(Acct.balance_ct);  // current balance
    Acct.m.PrintInDec("m"); 
    PrintSplitLine('-'); 
}

void SaveSP(SP &sp, std::string APGC_SP_File)
{
    std::ofstream fout;
    fout.open(APGC_SP_File, std::ios::binary); 
    fout << sp.ska;
    fout.close();   
}

void FetchSP(SP &sp, std::string APGC_SP_File)
{
    std::ifstream fin; 
    fin.open(APGC_SP_File, std::ios::binary); 
    fin >> sp.ska; 
    fin.close();   
}

void SavePP(PP &pp, std::string APGC_PP_File)
{
    std::ofstream fout; 
    fout.open(APGC_PP_File, std::ios::binary); 

    fout << pp.MAX_RECEIVER_NUM; 
    fout << pp.SN_LEN;
    fout << pp.MAXIMUM_COINS;  
    fout << pp.pka; 

    fout << pp.bullet_part; 
    fout << pp.enc_part; 

    fout.close();   
}

void FetchPP(PP &pp, std::string APGC_PP_File)
{
    std::ifstream fin; 
    fin.open(APGC_PP_File, std::ios::binary); 

    fin >> pp.MAX_RECEIVER_NUM;
    fin >> pp.SN_LEN; 
    fin >> pp.MAXIMUM_COINS;  
    fin >> pp.pka; 
 
    fin >> pp.bullet_part;
    fin >> pp.enc_part; 

    fin.close();   
}

void SaveAccount(Account &user, std::string APGC_Account_File)
{
    std::ofstream fout; 
    fout.open(APGC_Account_File, std::ios::binary);
    fout << user.identity;  
    fout << user.pk;              
    fout << user.sk;   
    fout << user.balance_ct;  
    fout << user.m; 
    fout.close();  
}

void FetchAccount(Account &user, std::string APGC_Account_File)
{
    std::ifstream fin; 
    fin.open(APGC_Account_File, std::ios::binary);
    fin >> user.identity; 
    fin >> user.pk;              
    fin >> user.sk;             
    fin >> user.balance_ct;
    fin >> user.m; 
    fin.close();  
}

/* This function implements Setup algorithm of APGC */
std::tuple<PP, SP> Setup(size_t LOG_MAXIMUM_COINS, size_t MAX_RECEIVER_NUM, size_t SN_LEN)
{
    PP pp; 
    SP sp; 

    pp.MAX_RECEIVER_NUM = MAX_RECEIVER_NUM; 
    if(IsPowerOfTwo(MAX_RECEIVER_NUM+1) == false){
        std::cerr << "parameters wrong: (MAX_RECEIVER_NUM+1) must be a power of 2" << std::endl; 
    }
    pp.SN_LEN = SN_LEN;    
    pp.MAXIMUM_COINS = BigInt(uint64_t(pow(2, LOG_MAXIMUM_COINS)));  


    size_t MAX_AGG_NUM = pp.MAX_RECEIVER_NUM + 1; 

    pp.bullet_part = Bullet::Setup(LOG_MAXIMUM_COINS, MAX_AGG_NUM); 
    
    size_t TRADEOFF_NUM = 7;
    pp.enc_part = TwistedExponentialElGamal::Setup(LOG_MAXIMUM_COINS, TRADEOFF_NUM);  

    std::tie(pp.pka, sp.ska) = TwistedExponentialElGamal::KeyGen(pp.enc_part);

    return {pp, sp};
}

/* initialize the encryption part for faster decryption */
void Initialize(PP &pp)
{
    std::cout << "initialize APGC >>>" << std::endl;  
    TwistedExponentialElGamal::Initialize(pp.enc_part); 
    PrintSplitLine('-'); 
}

/* create an account for input identity */
Account CreateAccount(PP &pp, std::string identity, BigInt &init_balance)
{
    Account newAcct;
    newAcct.identity = identity;

    std::tie(newAcct.pk, newAcct.sk) = TwistedExponentialElGamal::KeyGen(pp.enc_part); // generate a keypair

    newAcct.m = init_balance; 

    // initialize account balance with 0 coins
    BigInt r = Hash::StringToBigInt(newAcct.identity); 
    newAcct.balance_ct = TwistedExponentialElGamal::Enc(pp.enc_part, newAcct.pk, init_balance, r);

    #ifdef DEMO
        std::cout << identity << "'s APGC account creation succeeds" << std::endl;
        newAcct.pk.Print("pk"); 
        std::cout << identity << "'s initial balance = "; 
        newAcct.m.PrintInDec(); 
        std::cout << std::endl;
        PrintSplitLine('-'); 
    #endif 

    return newAcct;
}

/* reveal the balance */ 
BigInt RevealBalance(PP &pp, Account &Acct)
{
    return TwistedExponentialElGamal::Dec(pp.enc_part, Acct.sk, Acct.balance_ct); 
}


/* 
** support one to many transactions 
*/

// define the structure for confidential transaction
struct ToManyCTx{
    BigInt sn;                        // serial number: uniquely defines a transaction
    size_t k;                         // number of receivers
    // memo information
    //TwistedExponentialElGamal::CT sender_balance_ct;        // the current balance of pk1 (not necessarily included)
    //ECPoint pks;      // sender = pks
    //TwistedExponentialElGamal::CT sender_transfer_ct;
    std::vector<ECPoint> vec_pk;  
    std::vector<TwistedExponentialElGamal::MRCT> vec_participant_transfer_ct;    // (X0 = pka^r, X1 = pkr^r, Y = g^r h^v)
    std::vector<TwistedExponentialElGamal::CT> vec_participant_balance_ct;  

    ECPoint vector_commitment_B_l0;
    TwistedExponentialElGamal::CT refresh_updated_ct;
    // validity proof
    WellFormProduct::Proof plaintext_wellformed_proof;
    MutiliPlaintextEquality::Proof superivisor_plaintext_wellformed_proof;
    SumZero::Proof plaintext_sumzero_proof;
    AnyOutOfMany::Proof slack_participant_proof;
    Solvent_Equal::Proof solvent_equal_proof;
    Bullet::Proof bullet_right_solvent_proof;
    //TwistedExponentialElGamal::CT refresh_sender_updated_balance_ct;  // fresh encryption of updated balance (randomness is known)
    //PlaintextKnowledge::Proof plaintext_knowledge_proof; // NIZKPoK for refresh ciphertext (X^*, Y^*)
    //DLOGEquality::Proof correct_refresh_proof;     // fresh updated balance is correct

    std::vector<PlaintextEquality::Proof> vec_plaintext_equality_proof;     // NIZKPoK for transfer ciphertext (X1, X2, Y)
    //Bullet::Proof bullet_right_solvent_proof;   // aggregated range proof for v, m-v and v_i lie in the right range 

    DLOGKnowledge::Proof balance_proof; // prove v = v_1 +...+ v_n    
};

struct AnonSet{
    std::string identity;
    ECPoint pk;
    //size_t type; type = 0, sender, type = 1, receiver, type = 2: slack participant; 
    TwistedExponentialElGamal::CT balance_ct; // current balance
};

// save CTx into sn.ctx file
void SaveCTx(ToManyCTx &newCTx, std::string APGC_CTx_File)
{
    std::ofstream fout; 
    fout.open(APGC_CTx_File, std::ios::binary); 
    
    // save sn
    fout << newCTx.sn; 
     
    // save memo info
  
    for(auto i = 0; i < newCTx.vec_pk.size(); i++){
        fout << newCTx.vec_pk[i];
    }
    for(auto i = 0; i < newCTx.vec_participant_transfer_ct.size(); i++){
        fout << newCTx.vec_participant_transfer_ct[i];
    } 
    
    // save proofs
    // for(auto i = 0; i < newCTx.vec_plaintext_equality_proof.size(); i++){
    //     fout << newCTx.vec_plaintext_equality_proof[i];
    // }
    // fout << newCTx.refresh_sender_updated_balance_ct; 
    // fout << newCTx.correct_refresh_proof; 
    // fout << newCTx.bullet_right_solvent_proof; 
    // fout << newCTx.plaintext_knowledge_proof; 

    fout.close();

    // calculate the size of ctx_file
    std::ifstream fin; 
    fin.open(APGC_CTx_File, std::ios::ate | std::ios::binary);
    std::cout << APGC_CTx_File << " size = " << fin.tellg() << " bytes" << std::endl;
    fin.close(); 
}

std::string ExtractToSignMessageFromCTx(ToManyCTx &newCTx)
{
    // std::string str;
    // str += newCTx.sn.ToHexString() + newCTx.pks.ToByteString(); 
    // for(auto i = 0; i < newCTx.vec_pkr.size(); i++){
    //     str += newCTx.vec_pkr[i].ToByteString();
    // }

    // str += TwistedExponentialElGamal::CTToByteString(newCTx.sender_balance_ct);  
    // str += TwistedExponentialElGamal::CTToByteString(newCTx.sender_transfer_ct);  
    // for(auto i = 0; i < newCTx.vec_receiver_transfer_ct.size(); i++){
    //     str += TwistedExponentialElGamal::MRCTToByteString(newCTx.vec_receiver_transfer_ct[i]);
    // }

    // for(auto i = 0; i < newCTx.vec_plaintext_equality_proof.size(); i++){
    //     str += PlaintextEquality::ProofToByteString(newCTx.vec_plaintext_equality_proof[i]);  
    // }
   
    // str += Bullet::ProofToByteString(newCTx.bullet_right_solvent_proof);   
    // str += TwistedExponentialElGamal::CTToByteString(newCTx.refresh_sender_updated_balance_ct);  
    // str += PlaintextKnowledge::ProofToByteString(newCTx.plaintext_knowledge_proof); 

    // str += DLOGKnowledge::ProofToByteString(newCTx.balance_proof); 

    // return str;
}

/* 
* generate a confidential transaction: pks transfers vi coins to pkr[i] 
*/

size_t GetTheNthBit(size_t index, size_t n)
{
    return (index>>n)&1;
}

ToManyCTx CreateCTx(PP &pp, Account &Acct_sender, std::vector<BigInt> &vec_v, std::vector<ECPoint> &vec_pkr, std::vector<AnonSet> &vec_AnonSet, size_t sender_index, std::vector<size_t> vec_index)
{ 
    ToManyCTx newCTx; 
    newCTx.sn = GenRandomBigIntLessThan(order);
    int k = vec_v.size();
    size_t n = vec_AnonSet.size(); 
    if(k >= n){
        std::cerr << "the number of receivers must be less than total number" << n << std::endl;
    }
    if(k != vec_pkr.size()){
        std::cerr << "the number of receivers must be equal to the number of receiver's public key" << std::endl;
    }

    if(IsPowerOfTwo(n) == false){
        std::cerr << "receiver num must be 2^n-1" << std::endl;
    }

    std::string ctx_type = "(1-to-"+std::to_string(k)+")"; 

    #ifdef DEMO
        std::cout << "begin to genetate " << ctx_type << " ctx >>>>>>" << std::endl; 
    #endif
    PrintSplitLine('-'); 

    #ifdef DEMO
        std::cout <<"1. generate memo info of ctx" << std::endl;  
    #endif

    auto start_time = std::chrono::steady_clock::now();

    std::string transcript_str = "";  
    //newCTx.sn = Acct_sender.sn;
    newCTx.k = k;
    
    BigInt v = bn_0;
    for(auto i = 0; i < k; i++){
        v += vec_v[i]; 
    }
    // we will use the random position to generate the memo later
    // size_t sender_index;
    // std::vector<size_t> vec_index(k);
    // std::set<ECPoint> set_pkr;
    // for(auto i = 0; i < k; i++){
    //     set_pkr.insert(vec_pkr[i]);
    // }
    // for(auto i = 0; i < n; i++){
    //     if(vec_AnonSet[i].identity == Acct_sender.identity){
    //         sender_index = i;
    //     }
    //     else
    //     {
    //         if(set_pkr.find(vec_AnonSet[i].pk) != set_pkr.end()){
    //             vec_index.push_back(i);
    //         }
    //     }
    // }
    // std::sort(vec_index.begin(), vec_index.end());

    std::vector<ECPoint> vec_pk(n);
    std::vector<TwistedExponentialElGamal::MRCT> vec_participant_transfer_ct(n);
    std::vector<BigInt> vec_r(n);

    newCTx.vec_pk.resize(n);
    newCTx.vec_participant_transfer_ct.resize(n);
    std::vector<ECPoint> vec_pk_multi(2); 
    size_t cnt=0;
    for(auto i = 0; i < n; i++)
    {

       newCTx.vec_participant_balance_ct[i] = vec_AnonSet[i].balance_ct;
       if(i == sender_index)
       {
            newCTx.vec_pk[i] = Acct_sender.pk;
            vec_pk_multi[0] = Acct_sender.pk;
            vec_pk_multi[1] = pp.pka;
            vec_r[i] = GenRandomBigIntLessThan(order);
            newCTx.vec_participant_transfer_ct[0] = TwistedExponentialElGamal::Enc(pp.enc_part, vec_pk_multi, -v, vec_r[0]);
       }
       else
       {
            if(i == vec_index[cnt])
            {
                newCTx.vec_pk[i] = vec_AnonSet[i].pk;
                vec_pk_multi[0] = vec_AnonSet[i].pk;
                vec_pk_multi[1] = pp.pka;
                vec_r[i] = GenRandomBigIntLessThan(order);
                newCTx.vec_participant_transfer_ct[i] = TwistedExponentialElGamal::Enc(pp.enc_part, vec_pk_multi, vec_v[cnt], vec_r[i]);
                cnt++;
            }
            else
            {
                newCTx.vec_pk[i] = vec_AnonSet[i].pk;
                vec_pk_multi[0] = vec_AnonSet[i].pk;
                vec_pk_multi[1] = pp.pka;
                vec_r[i] = GenRandomBigIntLessThan(order);
                newCTx.vec_participant_transfer_ct[i] = TwistedExponentialElGamal::Enc(pp.enc_part, vec_pk_multi, bn_0, vec_r[i]);
            }
        }
    }
        
 
    #ifdef DEMO
        std::cout << "2. generate NIZKPoK for tx" << std::endl;  
    #endif
    WellFormProduct::PP plaintext_wellform_product_pp = WellFormProduct::Setup(pp.enc_part.g, pp.enc_part.h, n);
    WellFormProduct::Instance plaintext_wellform_product_instance;

    plaintext_wellform_product_instance.vec_pk = newCTx.vec_pk;
    plaintext_wellform_product_instance.vec_CL.resize(n);
    plaintext_wellform_product_instance.vec_CR.resize(n);
    for(auto i = 0; i < n; i++){
        plaintext_wellform_product_instance.vec_CL[i] = newCTx.vec_participant_transfer_ct[i].vec_X[0];
        plaintext_wellform_product_instance.vec_CR[i] = newCTx.vec_participant_transfer_ct[i].Y;
    }
    
    WellFormProduct::Witness plaintext_wellform_product_witness;

    plaintext_wellform_product_witness.vec_r = vec_r;
    std::vector<BigInt> vec_v_plaintext_wellform(n);
    cnt=0;
    for(auto i = 0; i < n; i++)
    {
        if(i == sender_index)
        {
            vec_v_plaintext_wellform[i] = -v;
        }
        else
        {
            if(i == vec_index[cnt])
            {
                vec_v_plaintext_wellform[i] = vec_v[cnt];
                cnt++;
            }
            else
            {
            vec_v_plaintext_wellform[i] = bn_0;
            }
        }
        
    }
    plaintext_wellform_product_witness.vec_v = vec_v_plaintext_wellform;

    std::string transcript_str = "";
    WellFormProduct::Proof plaintext_wellform_product_proof;
    WellFormProduct::Prove(plaintext_wellform_product_pp, plaintext_wellform_product_instance, 
                        plaintext_wellform_product_witness, transcript_str, plaintext_wellform_product_proof);

    SumZero::PP plaintext_sumzero_pp = SumZero::Setup(pp.enc_part.g, n);
    SumZero::Instance plaintext_sumzero_instance;
    plaintext_sumzero_instance.C = newCTx.vec_participant_transfer_ct[0].Y;
    for(auto i = 1; i < n; i++){
        plaintext_sumzero_instance.C += newCTx.vec_participant_transfer_ct[i].Y;
    }
  
    SumZero::Witness plaintext_sumzero_witness;
    plaintext_sumzero_witness.r = vec_r;
    transcript_str = "";
    SumZero::Proof plaintext_sumzero_proof;
    SumZero::Prove(plaintext_sumzero_pp, plaintext_sumzero_instance, plaintext_sumzero_witness, transcript_str, plaintext_sumzero_proof);


    #ifdef DEMO
        PrintSplitLine('-'); 
    #endif

    #ifdef DEMO
        std::cout << "2. generate NIZKPoK for slack_participant" << std::endl;  
    #endif
    
    AnyOutOfMany::PP slack_participant_pp = AnyOutOfMany::Setup(n, pp.enc_part.g, pp.enc_part.h);
    AnyOutOfMany::Instance slack_participant_instance;
    slack_participant_instance.vec_com.resize(n);
    for(auto i = 0; i < n; i++)
    {
        slack_participant_instance.vec_com[i] = newCTx.vec_participant_transfer_ct[i].Y;
    }
    AnyOutOfMany::Witness slack_participant_witness;
    //slack_participant_witness.vec_s.resize(pp.SLACK_PARTICIPANT_NUM);
    slack_participant_witness.vec_b.resize(n);
    for(auto i = 0; i < n; i++)
    {
        if(i != sender_index && std::find(vec_index.begin(), vec_index.end(), i) == vec_index.end())
        {
            slack_participant_witness.vec_b[i] = bn_1;
            slack_participant_witness.vec_s.push_back(vec_r[i]);
        }
        else
        {
            slack_participant_witness.vec_b[i] = bn_0;
        }  
    }
    transcript_str = "";
    AnyOutOfMany::Proof slack_participant_proof;
    AnyOutOfMany::Prove(slack_participant_pp, slack_participant_instance, slack_participant_witness, slack_participant_proof, transcript_str);
    
    #ifdef DEMO
        std::cout << "3. generate NIZKPoK for solvent" << std::endl;  
    #endif

    size_t m = log(n);
    std::vector<ECPoint> base_g = GenRandomECPointVector(m);
    ECPoint base_h = pp.enc_part.h;

    std::vector<BigInt> vec_l0;
    for(auto i = m-1; i >= 0; i--){
        if(GetTheNthBit(sender_index, i) == 1){
            vec_l0.push_back(bn_1);
        }
        else{
            vec_l0.push_back(bn_0);
        }
    }
    BigInt rb_l0 = GenRandomBigIntLessThan(order);
    ECPoint B_l0 = ECPointVectorMul(base_g, vec_l0) + base_h * rb_l0;
    newCTx.vector_commitment_B_l0 = B_l0;
    BigInt balance_sender = TwistedExponentialElGamal::Dec(pp.enc_part, Acct_sender.sk, Acct_sender.balance_ct);
    balance_sender = balance_sender - v;

    BigInt r_refresh = GenRandomBigIntLessThan(order);
    newCTx.refresh_updated_ct.X = Acct_sender.pk * r_refresh;
    newCTx.refresh_updated_ct.Y = pp.enc_part.g * r_refresh + pp.enc_part.h * balance_sender;

    Bullet::Instance bullet_instance_solvent;
    bullet_instance_solvent.C = {newCTx.refresh_updated_ct.Y};

    Bullet::Witness bullet_witness_solvent;  
    bullet_witness_solvent.r = {r_refresh}; 
    bullet_witness_solvent.v = {balance_sender};

    Bullet::Prove(pp.bullet_part, bullet_instance_solvent, bullet_witness_solvent, transcript_str, newCTx.bullet_right_solvent_proof); 

    std::vector<ECPoint> sum_ct_left(n);
    std::vector<ECPoint> sum_ct_right(n);
    for(auto i = 0; i < n; i++){
        sum_ct_left[i] = newCTx.vec_participant_transfer_ct[i].vec_X[0] + newCTx.vec_participant_balance_ct[i].X;
        sum_ct_right[i] = newCTx.vec_participant_transfer_ct[i].Y + newCTx.vec_participant_balance_ct[i].Y;
    }

    Solvent_Equal::PP solvent_equal_pp = Solvent_Equal::Setup(pp.enc_part.g, pp.enc_part.h, n);
    Solvent_Equal::Instance solvent_equal_instance;
    solvent_equal_instance.B = newCTx.vector_commitment_B_l0;
    solvent_equal_instance.Sum_CL = sum_ct_left;
    solvent_equal_instance.Sum_CR = sum_ct_right;
    solvent_equal_instance.Refresh_CL = newCTx.refresh_updated_ct.X;
    solvent_equal_instance.Refresh_CR = newCTx.refresh_updated_ct.Y;
    solvent_equal_instance.pk = newCTx.vec_pk;

    Solvent_Equal::Witness solvent_equal_witness;
    solvent_equal_witness.sk = Acct_sender.sk;
    solvent_equal_witness.l0 = sender_index;
    solvent_equal_witness.rb_l0 = rb_l0;
    solvent_equal_witness.r_refresh = r_refresh;
    solvent_equal_witness.balance_sender = balance_sender;

    Solvent_Equal::Proof solvent_equal_proof;

    transcript_str = "";

    Solvent_Equal::Prove(solvent_equal_pp, solvent_equal_instance, solvent_equal_witness, transcript_str, solvent_equal_proof);

    #ifdef DEMO
        std::cout << "4. generate NIZKPoK for supervise" << std::endl;  
    #endif

    MutiliPlaintextEquality::PP superivisor_plaintext_wellformed_pp = MutiliPlaintextEquality::Setup(pp.enc_part.g, pp.enc_part.h, n);
    MutiliPlaintextEquality::Instance superivisor_plaintext_wellformed_instance;

    superivisor_plaintext_wellformed_instance.pk_a = pp.pka;
    superivisor_plaintext_wellformed_instance.vec_CL.resize(n);
    superivisor_plaintext_wellformed_instance.vec_CR.resize(n);
    for(auto i = 0; i < n; i++){
        superivisor_plaintext_wellformed_instance.vec_CL[i] = newCTx.vec_participant_transfer_ct[i].vec_X[1];
        superivisor_plaintext_wellformed_instance.vec_CR[i] = newCTx.vec_participant_transfer_ct[i].Y;
    }

    MutiliPlaintextEquality::Witness superivisor_plaintext_wellformed_witness;

    superivisor_plaintext_wellformed_witness.vec_r = vec_r; // r is as same as the plaintext_wellform_product_witness
    superivisor_plaintext_wellformed_witness.vec_v = vec_v_plaintext_wellform;// v is as same as the plaintext_wellform_product_witness
    MutiliPlaintextEquality::Proof superivisor_plaintext_wellformed_proof;
    transcript_str = "";
    MutiliPlaintextEquality::Prove(superivisor_plaintext_wellformed_pp, superivisor_plaintext_wellformed_instance, 
                        superivisor_plaintext_wellformed_witness, transcript_str, superivisor_plaintext_wellformed_proof);

    auto end_time = std::chrono::steady_clock::now(); 

    auto running_time = end_time - start_time;
    std::cout << ctx_type << " ctx generation takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    return newCTx; 
}

/* check if the given confidential transaction is valid */ 
bool VerifyCTx(PP &pp, ToManyCTx &newCTx)
{     
    size_t n = newCTx.vec_pk.size(); 
    if(IsPowerOfTwo(n) == false){
        std::cerr << "receiver num must be 2^n-1" << std::endl;
    }
    std::string ctx_type = "(1-to-"+std::to_string(newCTx.k)+")"; 
    
    #ifdef DEMO
        std::cout << "begin to verify " <<ctx_type << " ctx >>>>>>" << std::endl; 
    #endif

    std::string transcript_str = "";

    auto start_time = std::chrono::steady_clock::now(); 

    // generate NIZK proof for validity of transfer              
    bool condition1 = true;      
    
    #ifdef DEMO
        if (condition1) std::cout << "NIZKPoK for plaintext equality accepts" << std::endl; 
        else std::cout << "NIZKPoK for plaintext equality rejects" << std::endl; 
    #endif

    // check V2
    bool condition2; 
    

    #ifdef DEMO
        if (condition2) std::cout << "NIZKPoK for refresh updated balance accepts" << std::endl; 
        else std::cout << "NIZKPoK for refresh updated balance rejects" << std::endl; 
    #endif

    // check range proof
    bool condition3; 

    

    #ifdef DEMO
        if (condition3) std::cout << "range proofs for transfer amount and updated balance accept" << std::endl; 
        else std::cout << "range proofs for transfer amount and updated balance reject" << std::endl; 
    #endif

    // check balance proof
    bool condition4;

    DLOGKnowledge::PP dlog_knowledge_pp = DLOGKnowledge::Setup();

    #ifdef DEMO
        if (condition4) std::cout << "NIZKPoK for balance proof accepts" << std::endl; 
        else std::cout << "NIZKPoK for balance proof rejects" << std::endl; 
    #endif

    // check the NIZK proof for refresh correctness
    bool condition5;

    #ifdef DEMO
        if (condition5) std::cout << "NIZKPoK for refreshing correctness accepts and ctx is authenticated" << std::endl; 
        else std::cout << "NIZKPoK for refreshing correctness rejects or ctx is unauthenticated" << std::endl; 
    #endif

    bool Validity = condition1 && condition2 && condition3 && condition4 && condition5; 

    std::string ctx_file = GetCTxFileName(newCTx); 
    #ifdef DEMO
        if (Validity) std::cout << ctx_file << " is valid <<<<<<" << std::endl; 
        else std::cout << ctx_file << " is invalid <<<<<<" << std::endl;
    #endif

    auto end_time = std::chrono::steady_clock::now(); 

    auto running_time = end_time - start_time;
    std::cout << ctx_type << " ctx verification takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    return Validity; 
}


/* print the details of a confidential to-many-transaction */
void PrintCTx(ToManyCTx &newCTx)
{
    // PrintSplitLine('-');
    // std::string ctx_file = GetCTxFileName(newCTx);  
    // std::cout << ctx_file << " content >>>>>>" << std::endl; 

    // std::cout << "current sender balance >>>" << std::endl; 
    // TwistedExponentialElGamal::PrintCT(newCTx.sender_balance_ct);
    // std::cout << std::endl; 

    // newCTx.pks.Print("sender's public key"); 
    // std::cout << "receiver's public key" << std::endl;
    // PrintECPointVector(newCTx.vec_pkr, "pkr");
    // std::cout << std::endl;  


    // std::cout << "sender's transfer ct >>>" << std::endl;
    // TwistedExponentialElGamal::PrintCT(newCTx.sender_transfer_ct);
    // std::cout << std::endl;

    // std::cout << "receiver's transfer ct >>>" << std::endl; 
    // for(auto i = 0; i < newCTx.vec_receiver_transfer_ct.size(); i++){
    //     TwistedExponentialElGamal::PrintCT(newCTx.vec_receiver_transfer_ct[i]); 
    //     std::cout << std::endl;
    // } 

    // std::cout << "NIZKPoK for plaintext equality of receiver's ct >>>" << std::endl; 
    // for(auto i = 0; i < newCTx.vec_plaintext_equality_proof.size(); i++){
    //     PlaintextEquality::PrintProof(newCTx.vec_plaintext_equality_proof[i]); 
    // }
    // std::cout << std::endl; 

    // std::cout << "NIZKPoK for input-output balance >>> " << std::endl; 
    // DLOGKnowledge::PrintProof(newCTx.balance_proof); // prove v = v_1 +...+ v_n 
    // std::cout << std::endl;

    // std::cout << "range proofs for transfer amount and updated balance >>> " << std::endl; 
    // Bullet::PrintProof(newCTx.bullet_right_solvent_proof); 
    // std::cout << std::endl;

    // std::cout << "refresh updated balance >>>" << std::endl;
    // TwistedExponentialElGamal::PrintCT(newCTx.refresh_sender_updated_balance_ct); 
    // std::cout << std::endl;

    // std::cout << "NIZKPoK of refresh updated balance >>>" << std::endl; 
    // PlaintextKnowledge::PrintProof(newCTx.plaintext_knowledge_proof); 
    // std::cout << std::endl;
    
    // std::cout << "NIZKPoK for refreshing correctness >>>" << std::endl; 
    // DLOGEquality::PrintProof(newCTx.correct_refresh_proof);     // fresh updated balance is correct
    // std::cout << std::endl;

    PrintSplitLine('-'); 
}

/* update Account if CTx is valid */
bool UpdateAccount(PP &pp, ToManyCTx &newCTx, Account &Acct_sender, std::vector<Account> &vec_Acct_participant)
{    
    //Acct_sender.sn = Acct_sender.sn + bn_1;

    // update sender's balance
    // Acct_sender.balance_ct = TwistedExponentialElGamal::HomoSub(Acct_sender.balance_ct, newCTx.sender_transfer_ct); 
    // Acct_sender.m = TwistedExponentialElGamal::Dec(pp.enc_part, Acct_sender.sk, Acct_sender.balance_ct); 
    // SaveAccount(Acct_sender, Acct_sender.identity+".account"); 

    TwistedExponentialElGamal::CT c_in; 
    for(auto i = 0; i < vec_Acct_participant.size(); i++){
        c_in.X = newCTx.vec_participant_transfer_ct[i].vec_X[0]; 
        c_in.Y = newCTx.vec_participant_transfer_ct[i].Y;
        // update participant's balance
        vec_Acct_participant[i].balance_ct = TwistedExponentialElGamal::HomoAdd(vec_Acct_participant[i].balance_ct, c_in); 
        vec_Acct_participant[i].m = TwistedExponentialElGamal::Dec(pp.enc_part, vec_Acct_participant[i].sk, vec_Acct_participant[i].balance_ct);
        SaveAccount(vec_Acct_participant[i], vec_Acct_participant[i].identity+".account"); 
    }

    return true; 
} 


/* check if a ctx is valid and update accounts if so */
bool Miner(PP &pp, ToManyCTx &newCTx, Account &Acct_sender, std::vector<Account> &vec_Acct_participant)
{
    for(auto i = 0; i < vec_Acct_participant.size(); i++){
        if (newCTx.vec_pk[i] != vec_Acct_participant[i].pk){
            std::cerr << i<<"-th participant does not match ctx" << std::endl; 
            return false;
        } 
    }

    std::string ctx_file = GetCTxFileName(newCTx); 
    if(VerifyCTx(pp, newCTx) == true){
        UpdateAccount(pp, newCTx, Acct_sender, vec_Acct_participant);
        SaveCTx(newCTx, ctx_file);  
        std::cout << ctx_file << " is recorded on the blockchain" << std::endl; 
        return true; 
    }
    else{
        std::cout << ctx_file << " is discarded" << std::endl; 
        return false; 
    }
}


/* supervisor opens CTx */
std::vector<BigInt> SuperviseCTx(SP &sp, PP &pp, ToManyCTx &ctx, std::vector<Account> &vec_Acct_participant)
{
    size_t n = ctx.vec_pk.size();
    std::vector<BigInt> vec_v(n); 

    std::cout << "Supervise " << GetCTxFileName(ctx) << std::endl; 
    auto start_time = std::chrono::steady_clock::now(); 

    //std::cout << ctx.pks.ToHexString() << " transfers "; 
    std::cout << std::endl;

    SuperviseResult result;
    TwistedExponentialElGamal::CT ct; 
    for(auto i = 0; i < n; i++){
        ct.X = ctx.vec_participant_transfer_ct[i].vec_X[1];
        ct.Y = ctx.vec_participant_transfer_ct[i].Y;  
        vec_v[i] = TwistedExponentialElGamal::Dec(pp.enc_part, sp.ska, ct);
        if(vec_v[i] < bn_0)
        {
            result.sender_pk = ctx.vec_pk[i];
        }
        else if(vec_v[i] > bn_0)
        {
            result.receiver_pks.push_back(ctx.vec_pk[i]);
            result.receiver_coins_values.push_back(vec_v[i]);
        }
        //std::cout << BN_bn2dec(vec_v[i].bn_ptr) << " coins to " << ctx.vec_pkr[i].ToHexString() << std::endl; 
    } 

    auto end_time = std::chrono::steady_clock::now(); 
    auto running_time = end_time - start_time;
    std::cout << "supervising ctx takes time = " 
    << std::chrono::duration <double, std::milli> (running_time).count() << " ms" << std::endl;

    return vec_v; 
}

}

#endif
