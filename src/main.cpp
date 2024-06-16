#include <mcl/bls12_381.hpp>
#include <iostream>
using namespace std;
using namespace mcl::bn;

#include "polynomial.hpp"
#include "proof.hpp"
#include "timer.hpp"

int main(int argc, char** argv) {
    initPairing(mcl::BLS12_381);

    uint deg = stoi(argv[1]);
    const auto pp = trustedSetup(deg);
    const auto& gVec = pp.gVec;
    const auto& hVec = pp.hVec;
    const auto& g2 = pp.g2;
    const auto& g2_tau = pp.g2_tau;

    Polynomial<Fr> F = randomPolynomial(deg);
    auto com_F = commitPoly(F, gVec);

    Fr x, rx, y, ry;
    x.setByCSPRNG();
    rx.setByCSPRNG();
    y = F(x);
    ry.setByCSPRNG();
    const auto& g = gVec[0];

    const auto& h = hVec[0];
    auto com_x = g * x + h * rx;
    auto com_y = g * y + h * ry;

    Timer ptimer, vtimer;
    bool accepted = SecretEval(y, ry, x, rx, F, com_y, com_x, com_F, pp, ptimer, vtimer);

    cout << "Accepted: " << accepted << endl;
    cout << "Proving time: " << ptimer.getTotalTime() << " s" << endl;
    cout << "Verification time: " << vtimer.getTotalTime() << " s" << endl;

    // Fr temp = getRootOfUnity(8);
    // Fr temp2;
    // Fr::pow(temp2, temp, 255);
    // cout << temp2 << endl;
    // Fr::pow(temp2, temp, 256);
    // cout << temp2 << endl;

    return 0;
}