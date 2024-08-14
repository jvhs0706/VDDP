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
    cout << len << endl;
    auto pp = LegendrePRNGTrustedSetup(len);

    vector<Fr> va(len), vb(len), vc(len), vd(len);
    vector<bool> vs(len);

    for (uint i = 0; i < len; ++ i)
    {
        va[i].setByCSPRNG();
        vb[i].setByCSPRNG();
        vc[i] = va[i] * vb[i];
        vs[i] = rand() % 2;
        vd[i] = vs[i] ? va[i] : vc[i];
    }

    const auto& omega = pp.omega_gen;


    Polynomial Fa = ntt_vec_to_poly_given_omega(va, omega);
    Polynomial Fb = ntt_vec_to_poly_given_omega(vb, omega);
    Polynomial Fc = ntt_vec_to_poly_given_omega(vc, omega);
    Polynomial Fd = ntt_vec_to_poly_given_omega(vd, omega);
    Polynomial Fs = ntt_vec_to_poly_given_omega(vector<Fr>(vs.begin(), vs.end()), omega);

    Polynomial Ra = randomPolynomial(Fa.getDegree());
    Polynomial Rb = randomPolynomial(Fb.getDegree());
    Polynomial Rc = randomPolynomial(Fc.getDegree());
    Polynomial Rd = randomPolynomial(Fd.getDegree());
    Polynomial Rs = randomPolynomial(Fs.getDegree());

    G1 com_a = commitPoly(Fa, Ra, pp.pp.gVec, pp.pp.hVec);
    G1 com_b = commitPoly(Fb, Rb, pp.pp.gVec, pp.pp.hVec);
    G1 com_c = commitPoly(Fc, Rc, pp.pp.gVec, pp.pp.hVec);
    G1 com_d = commitPoly(Fd, Rd, pp.pp.gVec, pp.pp.hVec);
    G1 com_s = commitPoly(Fs, Rs, pp.pp.gVec, pp.pp.hVec);

    Timer ptimer, vtimer;
    auto accepted = Hadamard(
        Fc, Rc, Fa, Ra, Fb, Rb, 
        len, omega, com_c, com_a, com_b, 
        pp.pp, ptimer, vtimer
    );
    cout << "Accepted: " << accepted << endl;

    cout << "Proving time: " << ptimer.getTotalTime() << endl;
    cout << "Verification time: " << vtimer.getTotalTime() << endl;

    Timer ptimer2, vtimer2;
    accepted = Mux(
        Fs, Rs, Fa, Ra, Fb, Rb, Fd, Rd, 
        len, omega, com_s, com_a, com_b, com_d, 
        pp.pp, ptimer2, vtimer2
    );
    cout << "Accepted: " << accepted << endl;
    cout << "Proving time: " << ptimer2.getTotalTime() << endl;
    cout << "Verification time: " << vtimer2.getTotalTime() << endl;
    





    // vector<Fr> rt_vec;

    // Fr key, r_key;
    // key.setByCSPRNG();
    // r_key.setByCSPRNG();
    
    // auto random_bits = LegendrePRNG(key, len, rt_vec);
    // Polynomial F_rt, F_res;
    // convertLegendrePRNG(rt_vec, random_bits, F_rt, F_res, pp);

    // Polynomial R_rt = randomPolynomial(len), R_res = randomPolynomial(len);
    // G1 com_key, com_rt, com_res;
    // commitLegendrePRNG(key, F_rt, F_res, r_key, R_rt, R_res, com_key, com_rt, com_res, pp);

    // Timer ptimer, vtimer;
    // auto accepted = proveLegendrePRNG(key, F_rt, F_res, 
    //     r_key, R_rt, R_res,
    //     com_key, com_rt, com_res, 
    //     pp, ptimer, vtimer);

    // cout << "Accepted: " << accepted << endl;

    return 0;
}