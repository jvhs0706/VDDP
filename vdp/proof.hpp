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

G1 commitPoly(const Polynomial& F, const vector<G1>& gVec);

G1 commitPoly(const Polynomial& F, const Polynomial& R, const vector<G1>& gVec, const vector<G1>& hVec);

G1 provePolyEval(const Polynomial& F, 
    const Fr& x, Fr& y, const vector<G1>& gVec);

G1 provePolyEval(const Polynomial& F, const Polynomial& R, 
    const Fr& x, 
    Fr& yF, Fr& yR,
    const vector<G1>& gVec, const vector<G1>& hVec);

bool verifyPolyEval(const Fr& yF, const Fr& yR, const G1& proof, const G1& commitment, const Fr& x, const G1& g, const G1& h, const G2& g2, const G2& g2_tau);

bool verifyPolyEval(const Fr& y, const G1& proof, const G1& commitment, const Fr& x, const G1& g, const G2& g2, const G2& g2_tau);

// Evaluating a public polynomial at a secret point
bool SecretEval(const Fr& y, const Fr& ry, const Fr& x, const Fr& rx, const Polynomial& F, 
    const G1& com_y, const G1& com_x, const G1& com_F, 
    const PubParam& pp, Timer& ptimer, Timer& vtimer);

// Evaluating a secret polynomial at a public point
bool EvalSecret(const Fr& y, const Fr& ry, const Fr& x, // secret y, public x
    const Polynomial& F, const Polynomial& R, // secret F, R
    const G1& com_y, const G1& com_F, 
    const PubParam& pp, Timer& ptimer, Timer& vtimer);

bool Equal(const Fr& x, const Fr& r_x, const Fr& y, const Fr& r_y, 
    const G1& com_x, const G1& com_y, const G1& g, const G1& h,
    Timer& ptimer, Timer& vtimer);

// z == x * y
bool Prod(const Fr& z, const Fr& r_z, const Fr& x, const Fr& r_x, const Fr& y, const Fr& r_y,
    const G1& com_z, const G1& com_x, const G1& com_y, const G1& g, const G1& h,
    Timer& ptimer, Timer& vtimer);

//The followings are for the batched versions

bool Binary(const Polynomial& F, const Polynomial& R, 
    const uint len, const Fr& omega, 
    const G1& com_F, const PubParam& pp, 
    Timer& ptimer, Timer& vtimer);

#endif // PROOF_HPP