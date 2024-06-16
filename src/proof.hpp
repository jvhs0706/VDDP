#ifndef PROOF_HPP
#define PROOF_HPP

#include <mcl/bls12_381.hpp>
#include <vector>
#include "polynomial.hpp"
#include "timer.hpp"

using namespace std;
using namespace mcl::bn;

struct PubParam {
    vector<G1> gVec;
    vector<G1> hVec;
    G2 g2;
    G2 g2_tau;
};

PubParam trustedSetup(uint maxDeg);

G1 commitPoly(const Polynomial<Fr>& F, const vector<G1>& gVec);

G1 commitPoly(const Polynomial<Fr>& F, const Polynomial<Fr>& R, const vector<G1>& gVec, const vector<G1>& hVec);

G1 provePolyEval(const Polynomial<Fr>& F, 
    const Fr& x, Fr& y, const vector<G1>& gVec);

G1 provePolyEval(const Polynomial<Fr>& F, const Polynomial<Fr>& R, 
    const Fr& x, 
    Fr& yF, Fr& yR,
    const vector<G1>& gVec, const vector<G1>& hVec);

bool verifyPolyEval(const Fr& yF, const Fr& yR, const G1& proof, const G1& commitment, const Fr& x, const G1& g, const G1& h, const G2& g2, const G2& g2_tau);

bool verifyPolyEval(const Fr& y, const G1& proof, const G1& commitment, const Fr& x, const G1& g, const G2& g2, const G2& g2_tau);

bool SecretEval(const Fr& y, const Fr& ry, const Fr& x, const Fr& rx, const Polynomial<Fr>& F, 
    const G1& com_y, const G1& com_x, const G1& com_F, 
    const PubParam& pp, Timer& ptimer, Timer& vtimer);

bool Prod(const Fr& z, const Fr& r_z, const Fr& x, const Fr& r_x, const Fr& y, const Fr& r_y,
    const G1& com_x, const G1& com_y, const G1& com_z, const G1& g, const G1& h,
    Timer& ptimer, Timer& vtimer);

Polynomial<Fr> randomPolynomial(uint deg);

Fr getRootOfUnity(uint n);

#endif // PROOF_HPP