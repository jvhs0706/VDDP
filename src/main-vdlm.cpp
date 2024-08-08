#include <mcl/bls12_381.hpp>
#include <iostream>
using namespace std;
using namespace mcl::bn;

#include "utils.hpp"
#include "polynomial.hpp"
#include "proof.hpp"
#include "timer.hpp"
#include "prng.hpp"
#include "vdlm.hpp"

int main(int argc, char** argv) {
    initPairing(mcl::BLS12_381);
    // initialize the random seed
    srand(time(NULL));

    uint len = stoi(argv[1]);
    auto pp = LegendrePRNGTrustedSetup(len);

    vector<Fr> rt_vec;

    Fr key;
    key.setByCSPRNG();
    
    auto random_bits = LegendrePRNG(key, len, rt_vec);
    Polynomial F_rt, F_res;
    convertLegendrePRNG(rt_vec, random_bits, F_rt, F_res, pp);

    Polynomial R_rt = randomPolynomial(len), R_res = randomPolynomial(len);
    G1 com_rt, com_res;
    commitLegendrePRNG(F_rt, F_res, R_rt, R_res, com_rt, com_res, pp);

    Timer ptimer, vtimer;
    Binary(F_res, R_res, len, pp.omega_gen, com_res, pp.pp, ptimer, vtimer); 

    return 0;
}