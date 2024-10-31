# APGC: Anonymous privacy-preserving transaction system based on Kunlun---A Modern Crypto Library


## Overview

I give a C++ wrapper for OpenSSL, making it handy to use, without worrying about the cumbersome memory management and memorizing the complex interfaces. Based on this wrapper, I am going to build an efficient and modular crypto library. 


## Specifications

- OS: MAC OS x64, Linux x64
- Language: C++
- Requires: OpenSSL, OpenMP

## Install Depedent Libaraies
### On MACOS
* download the latest OpenSSL from the website, to support curve25519, 
modify crypto/ec/curve25519.c line 211: remove "static", then compile it:
```
  $ ./Configure darwin64-x86_64-cc shared enable-ec_nistp_64_gcc_128 no-ssl2 no-ssl3 no-comp --openssldir=/usr/local/ssl/macos-x86_64
  $ make depend
  $ sudo make install
```

test if the function x25519_scalar_mulx is available
```
  $ cd /usr/local/lib
  $ nm libcrypto.a | grep x25519_scalar_mulx
```

* install OpenMP
```
  $ brew install libomp 
```

<!-- * install abseil-cpp
```
  $ git clone git@github.com:abseil/abseil-cpp.git 
  $ mkdir build && cd build
  $ cmake -DABSL_BUILD_TESTING=ON -DABSL_USE_GOOGLETEST_HEAD=ON -DCMAKE_CXX_STANDARD=14 ..
  $ make install
``` -->


### On Linux
* install OpenSSL 3.0

do the same modification as in MACOS, then compile it according to
```
  $ ./Configure no-shared enable-ec_nistp_64_gcc_128 no-ssl2 no-ssl3 no-comp --prefix=/usr/local/openssl
  $ make depend
  $ sudo make install
```

if reporting cannot find "opensslv.h" error, try to install libssl-dev
```
  $ sudo apt-get install libssl-dev 
```

* install OpenMP
```
  $ sudo apt-get install libomp-dev 
```

## Code Structure

- README.md

- CmakeLists.txt: cmake file

- /build

- /include
  * std.inc: standard header files
  * openssl.inc: openssl header files
  * global.hpp: define global variables for kunlun lib as well as error information

- /utility: dependent files
  * bit_operation.hpp
  * routines.hpp: related routine algorithms 
  * print.hpp: print info for debug
  * murmurhash3.hpp: add fast non-cryptographic hash
  * polymul.hpp: naive poly mul
  * serialization.hpp: overload serialization for uint and string type data

- /crypto: C++ wrapper for OpenSSL
  * setup.hpp: initialize crypto environments, including big number, elliptic curves, and aes
  * ec_group.hpp: initialize ec group environment, define compressed-point on-off, precomputation on-off 
  * ec_point.hpp: class for EC_POINT of ordinary EC curves 
  * ec_25519.hpp: class for x25519 method of specific Curve25519 
  * bigint.hpp: class for BIGNUM, also include initialization of big num
  * hash.hpp: all kinds of cryptographic hash functions
  * aes.hpp: implement AES using SSE, as well as initialization of aes
  * prg.hpp: implement PRG associated algorithms
  * prp.hpp: implement PRP using AES
  * block.hpp: __m128i related algorithms (necessary for exploiting SSE)

- /pke: public key encryption schemes
  * twisted_exponential_elgamal.hpp
  * elgamal.hpp: standard ElGamal PKE whose message space is G 
  * calculate_dlog.hpp: implement optimized general Shank's algorithm

- /commitment
  * pedersen.hpp: multi-element Pedersen commitment

- /cryptocurrency
  * apgc.hpp: the apgc system 

- zkp
  - /apgcproofs
    * apgc_amorhom_proof.hpp: apgc_sdr_solvent_equal.hpp invoke this
    * apgc_eqmdl_proof.hpp: prove G=g^a and H = h^a and P = p^a
    * apgc_eqmdl4two_proof.hpp: prove P=g^a and H =h^a
    * apgc_sdr_solvent_equal.hpp: prove sender's balance is solvent
    * apgc_sum_zero.hpp: the sum of the commitment is zero
    * apgc_wellform.hpp: cipher is well-form
  - /nizk: associated sigma protocol for twisted elgamal; obtained via Fiat-Shamir transform  
    * nizk_1oon_for_sdr_trans.hpp: nizk_sdr_trans.hpp invoke this
    * nizk_1oon.hpp: nizk_1oon_for_sdr_trans.hpp invoke this
    * nizk_any_out_of_many.hpp: ant out of n is encryption of 0
    * nizk_lin_bit.hpp: vector bit commitment
    * nizk_sdr_trans.hpp: prove one-out-of-n ciphertexts is encryption of -v

  - /bulletproofs
    * bullet_proof.hpp: the aggregating logarithmic size bulletproofs
    * innerproduct_proof.hpp: the inner product argument (used by Bulletproof to shrink the proof size) 


---

## Compile and Run
```
  $ mkdir build && cd build
  $ cmake ..
  $ make
  $ ./test_apgc
```

---

## Multi-threads Support
- Kunlun supports multithread by leveraging openmp. The underlying OpenSSL is not thread-safe, cause several threads may access a critial data structure "bn_ctx" concurrently. Kunlun is made thread-safe by introducing an array of bn_ctx. Thus, each thread has its own bn_ctx.     

- The global setting for multi-thread support lies at "include/global.hpp" line 19

- For multi-thread (n)
```
inline const size_t NUMBER_OF_THREADS = n; 
the default value of n is NUMBER_OF_PHYSICAL_CORES 
```

- For single-thread
```
inline const size_t NUMBER_OF_THREADS = 1; 
```

## Elliptic curve setting
- Kunlun supports all EC curves provided by OpenSSL. The global setting of EC curves lies at "crypto/ec_group.hpp" line 16-18. 

```
inline int curve_id = NID_X9_62_prime256v1; // choose other curves by specifying curve-ID  
#define ECPOINT_COMPRESSED                  // comment this line to enable uncompressed representation
#define ENABLE_X25519_ACCELERATION      // (un)comment this line to enable x25519 acceleration method
```

Note: x22519 is an efficnet DDH-based non-interactive key exchange (NIKE) protocol based on curve25519. The essense of x25519 is exactly cwPRF. Its remarkable efficency is attained by performing "somehow EC exponentiation" with only X-coordinates (perhaps x25519 name after it). However, in x25519 the EC exponetiation is not standard, and EC addition is not well-defined. We stress that curve25519 certainly support standard EC exponentiation and addition, but x25519 method does not. Kunlun provides the option of using x25519 method to improve performance of applications when it is applicable (involving only cwPRF). But, since x25519 method is not full-fledged, ordinary EC curves are always necessary for base Naor-Pinkas OT. Therefore, users must specify one ordinary EC curve when implementing ECC.    


## License

This library is licensed under the [MIT License](LICENSE).

---



