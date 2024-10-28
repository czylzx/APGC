#include "../apgc/apgc.hpp"
#include "../crypto/setup.hpp"
#include <time.h>
size_t getrandom(size_t n)
{
    srand(time(0));
    return rand() % n;
}
void Build_APGC_Test_Enviroment(size_t number)
{
    PrintSplitLine('-'); 
    std::cout << "build test enviroment for APGC >>>" << std::endl; 
    PrintSplitLine('-'); 
    std::cout << "setup APGC system" << std::endl; 
    // setup APGC system
    size_t SN_LEN = 4;
    size_t LOG_MAXIMUM_COINS = 32;      
    size_t MAX_RECEIVER_NUM = 7;  

    APGC::SP sp;
    APGC::PP pp;

    std::tie(pp, sp) = APGC::Setup(LOG_MAXIMUM_COINS, MAX_RECEIVER_NUM, SN_LEN); 

    APGC::Initialize(pp);

    std::string APGC_SP_Filename = "APGC.sp"; 
    APGC::SaveSP(sp, APGC_SP_Filename); 

    std::string APGC_PP_Filename = "APGC.pp"; 
    APGC::SavePP(pp, APGC_PP_Filename); 

    std::cout << "press any key to continue >>>" << std::endl; 
    system ("read");


    // create accounts for Alice and Bob and Tax
    std::cout << "generate four accounts" << std::endl; 
    PrintSplitLine('-'); 

    BigInt Alice_balance = BigInt(512); 
    APGC::Account Acct_Alice = APGC::CreateAccount(pp, "Alice", Alice_balance); 
    std::string Alice_Acct_FileName = "Alice.account"; 
    APGC::SaveAccount(Acct_Alice, Alice_Acct_FileName); 

    BigInt Bob_balance = BigInt(256); 
    APGC::Account Acct_Bob = APGC::CreateAccount(pp, "Bob", Bob_balance); 
    std::string Bob_Acct_FileName = "Bob.account"; 
    APGC::SaveAccount(Acct_Bob, Bob_Acct_FileName); 

    BigInt Carl_balance = BigInt(128); 
    APGC::Account Acct_Carl = APGC::CreateAccount(pp, "Carl", Carl_balance); 
    std::string Carl_Acct_FileName = "Carl.account"; 
    APGC::SaveAccount(Acct_Carl, Carl_Acct_FileName); 

    BigInt David_balance = BigInt(64);
    APGC::Account Acct_David = APGC::CreateAccount(pp, "David", David_balance);
    std::string David_Acct_FileName = "David.account";
    APGC::SaveAccount(Acct_David, David_Acct_FileName);

    BigInt Eve_balance = BigInt(32);
    APGC::Account Acct_Eve = APGC::CreateAccount(pp, "Eve", Eve_balance);
    std::string Eve_Acct_FileName = "Eve.account";
    APGC::SaveAccount(Acct_Eve, Eve_Acct_FileName);

    BigInt Frank_balance = BigInt(16);
    APGC::Account Acct_Frank = APGC::CreateAccount(pp, "Frank", Frank_balance);
    std::string Frank_Acct_FileName = "Frank.account";
    APGC::SaveAccount(Acct_Frank, Frank_Acct_FileName);

    BigInt Grace_balance = BigInt(8);   
    APGC::Account Acct_Grace = APGC::CreateAccount(pp, "Grace", Grace_balance);
    std::string Grace_Acct_FileName = "Grace.account";
    APGC::SaveAccount(Acct_Grace, Grace_Acct_FileName);

    BigInt Henry_balance = BigInt(4);
    APGC::Account Acct_Henry = APGC::CreateAccount(pp, "Henry", Henry_balance);
    std::string Henry_Acct_FileName = "Henry.account";
    APGC::SaveAccount(Acct_Henry, Henry_Acct_FileName);

    BigInt Ida_balance = BigInt(32);
    APGC::Account Acct_Ida = APGC::CreateAccount(pp, "Ida", Ida_balance);
    std::string Ida_Acct_FileName = "Ida.account";
    APGC::SaveAccount(Acct_Ida, Ida_Acct_FileName);

    BigInt Jack_balance = BigInt(32);
    BigInt Jack_sn = bn_1;
    APGC::Account Acct_Jack = APGC::CreateAccount(pp, "Jack", Jack_balance);
    std::string Jack_Acct_FileName = "Jack.account";
    APGC::SaveAccount(Acct_Jack, Jack_Acct_FileName);

    BigInt Kate_balance = BigInt(32);
    APGC::Account Acct_Kate = APGC::CreateAccount(pp, "Kate", Kate_balance);
    std::string Kate_Acct_FileName = "Kate.account";
    APGC::SaveAccount(Acct_Kate, Kate_Acct_FileName);

    BigInt Leo_balance = BigInt(32);
    APGC::Account Acct_Leo = APGC::CreateAccount(pp, "Leo", Leo_balance);
    std::string Leo_Acct_FileName = "Leo.account";
    APGC::SaveAccount(Acct_Leo, Leo_Acct_FileName);

    BigInt Mary_balance = BigInt(32);
    APGC::Account Acct_Mary = APGC::CreateAccount(pp, "Mary", Mary_balance);
    std::string Mary_Acct_FileName = "Mary.account";
    APGC::SaveAccount(Acct_Mary, Mary_Acct_FileName);

    BigInt Nick_balance = BigInt(32);
    APGC::Account Acct_Nick = APGC::CreateAccount(pp, "Nick", Nick_balance);
    std::string Nick_Acct_FileName = "Nick.account";
    APGC::SaveAccount(Acct_Nick, Nick_Acct_FileName);

    BigInt Olivia_balance = BigInt(32);
    APGC::Account Acct_Olivia = APGC::CreateAccount(pp, "Olivia", Olivia_balance);
    std::string Olivia_Acct_FileName = "Olivia.account";
    APGC::SaveAccount(Acct_Olivia, Olivia_Acct_FileName);

    BigInt Paul_balance = BigInt(32);
    APGC::Account Acct_Paul = APGC::CreateAccount(pp, "Paul", Paul_balance);
    std::string Paul_Acct_FileName = "Paul.account";
    APGC::SaveAccount(Acct_Paul, Paul_Acct_FileName);

    // BigInt Tax_balance = bn_0; 
    // APGC::Account Acct_Tax = APGC::CreateAccount(pp, "Tax", Tax_balance); 
    // std::string Tax_Acct_FileName = "Tax.account"; 
    // APGC::SaveAccount(Acct_Tax, Tax_Acct_FileName);

    std::cout << "press any key to continue >>>" << std::endl; 
    system ("read");
} 

void Emulate_APGC_System(size_t number, size_t kreceiver)
{
    size_t RANGE_LEN = 32; // set the range to be [0, 2^32-1]
    size_t AGG_NUM = 2; 
    
    APGC::SP sp;  
    APGC::FetchSP(sp, "APGC.sp"); 

    APGC::PP pp;  
    APGC::FetchPP(pp, "APGC.pp"); 
    APGC::PrintPP(pp); 

    APGC::Account Acct_Alice;  
    APGC::FetchAccount(Acct_Alice, "Alice.account"); 
    APGC::PrintAccount(Acct_Alice); 

    APGC::Account Acct_Bob;  
    APGC::FetchAccount(Acct_Bob, "Bob.account"); 
    APGC::PrintAccount(Acct_Bob); 

    APGC::Account Acct_Carl;  
    APGC::FetchAccount(Acct_Carl, "Carl.account"); 
    APGC::PrintAccount(Acct_Carl); 

    APGC::Account Acct_Tax;  
    APGC::FetchAccount(Acct_Tax, "Tax.account"); 
    APGC::PrintAccount(Acct_Tax); 

    APGC::Account Acct_David;
    APGC::FetchAccount(Acct_David, "David.account");
    APGC::PrintAccount(Acct_David);

    APGC::Account Acct_Eve;
    APGC::FetchAccount(Acct_Eve, "Eve.account");
    APGC::PrintAccount(Acct_Eve);

    APGC::Account Acct_Frank;
    APGC::FetchAccount(Acct_Frank, "Frank.account");
    APGC::PrintAccount(Acct_Frank);

    APGC::Account Acct_Grace;
    APGC::FetchAccount(Acct_Grace, "Grace.account");
    APGC::PrintAccount(Acct_Grace);

    APGC::Account Acct_Henry;
    APGC::FetchAccount(Acct_Henry, "Henry.account");
    APGC::PrintAccount(Acct_Henry);

    APGC::Account Acct_Ida;
    APGC::FetchAccount(Acct_Ida, "Ida.account");
    APGC::PrintAccount(Acct_Ida);

    APGC::Account Acct_Jack;
    APGC::FetchAccount(Acct_Jack, "Jack.account");
    APGC::PrintAccount(Acct_Jack);

    APGC::Account Acct_Kate;
    APGC::FetchAccount(Acct_Kate, "Kate.account");
    APGC::PrintAccount(Acct_Kate);

    APGC::Account Acct_Leo;
    APGC::FetchAccount(Acct_Leo, "Leo.account");
    APGC::PrintAccount(Acct_Leo);

    APGC::Account Acct_Mary;
    APGC::FetchAccount(Acct_Mary, "Mary.account");
    APGC::PrintAccount(Acct_Mary);

    APGC::Account Acct_Nick;
    APGC::FetchAccount(Acct_Nick, "Nick.account");
    APGC::PrintAccount(Acct_Nick);

    APGC::Account Acct_Olivia;
    APGC::FetchAccount(Acct_Olivia, "Olivia.account");
    APGC::PrintAccount(Acct_Olivia);

    APGC::Account Acct_Paul;
    APGC::FetchAccount(Acct_Paul, "Paul.account");
    APGC::PrintAccount(Acct_Paul);

    std::vector<APGC::Account> vec_Acct= {Acct_Alice, Acct_Bob, Acct_Carl, Acct_David, 
                                                    Acct_Eve, Acct_Frank, Acct_Grace, Acct_Henry, 
                                                    Acct_Ida, Acct_Jack, Acct_Kate, Acct_Leo, Acct_Mary, 
                                                    Acct_Nick, Acct_Olivia, Acct_Paul};

    std::cout << "Alice is going to transfer >>>>>>>>" << std::endl;
    
    std::cout << "begin to the test of 1-to-" << kreceiver << " ctx" << std::endl; 
    PrintSplitLine('-'); 
    std::vector<APGC::AnonSet> vec_AnonSet(number);
    std::vector<BigInt> vec_v(kreceiver); 
    for(size_t i = 0; i < kreceiver; i++)
    {
        size_t value = getrandom(32);
        vec_v[i] = BigInt(value);
    }
    
    // std::cout << "Tax is going to transfer "<< std::endl; 
    // std::cout << BN_bn2dec(vec_v[0].bn_ptr) << " coins to Alice" << std::endl; 
    // std::cout << BN_bn2dec(vec_v[1].bn_ptr) << " coins to Bob" << std::endl; 
    // std::cout << BN_bn2dec(vec_v[2].bn_ptr) << " coins to Carl" << std::endl;

    /*Todo*/
    // std::vector<size_t> vec_index(kreceiver);
    // size_t sender_index = getrandom(number);
    // for(size_t i = 0; i < kreceiver; i++)
    // {
    //     size_t receiver_index = getrandom(number);
    //     while(receiver_index == sender_index)
    //     {
    //         receiver_index = getrandom(number);
    //     }
    //     vec_index[i] = receiver_index;
    // }
    //std::sort(vec_index.begin(), vec_index.end());
   
    std::vector<size_t> vec_index(kreceiver);
    size_t sender_index = 0;
    for(size_t i = 0; i < kreceiver; i++)
    {
        size_t receiver_index = i+1;
        vec_index[i] = receiver_index;
    }
    std::vector<ECPoint> vec_pkr(kreceiver);
    std::vector<APGC::Account> vec_Acct_participant(kreceiver);
    for(size_t i = 0; i < kreceiver; i++)
    {
        vec_pkr[i] = vec_Acct[vec_index[i]].pk;
        
    }
    for(size_t i = 0; i < number; i++)
    {
        vec_AnonSet[i].identity = vec_Acct[i].identity;
        vec_AnonSet[i].pk = vec_Acct[i].pk;
        vec_AnonSet[i].balance_ct = vec_Acct[i].balance_ct;
        vec_Acct_participant[i] = vec_Acct[i];
    }
    APGC::ToManyCTx ctx = APGC::CreateCTx(pp, Acct_Alice, vec_v, vec_pkr,vec_AnonSet, sender_index, vec_index);
    APGC::Miner(pp, ctx, Acct_Alice, vec_Acct_participant); 
    PrintSplitLine('-'); 

    std::cout << "after 1st valid 1-to-n ctx >>>>>>" << std::endl; 
    // PrintSplitLine('-'); 
    // APGC::PrintAccount(Acct_Alice); 
    // APGC::PrintAccount(Acct_Bob); 
    // APGC::PrintAccount(Acct_Carl); 
    // APGC::PrintAccount(Acct_Tax); 

    std::cout << "press any key to continue >>>" << std::endl; 
    system ("read");

    std::cout << "supervision of 1-to-n ctx begins >>>" << std::endl; 
    PrintSplitLine('-'); 
    APGC::SuperviseCTx(sp, pp, ctx, vec_Acct_participant); 
    PrintSplitLine('-'); 
    std::cout << "supervision of 1-to-k ctx ends >>>" << std::endl; 
    PrintSplitLine('-');

}



int main()
{
    CRYPTO_Initialize();   
    size_t number = 8;
    size_t kreceiver = 3;

    Build_APGC_Test_Enviroment(number); 
    Emulate_APGC_System(number, kreceiver);

    CRYPTO_Finalize(); 

    return 0; 
}



