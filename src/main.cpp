#include <mcl/bls12_381.hpp>
#include <iostream>
using namespace std;
using namespace mcl::bn;

#include "polynomial.hpp"
#include "poly-commit.hpp"
#include "timer.hpp"

Polynomial<Fr> randomPolynomial(uint deg) {
    vector<Fr> coefs;
    for (int i = 0; i < deg + 1; i++) {
        Fr c;
        c.setByCSPRNG();
        coefs.push_back(c);
    }
    return coefs;
}

int main(int argc, char** argv) {
    initPairing(mcl::BLS12_381);

    G1 P;
    const char *g1Str = "1 0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb 0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1";
    P.setStr(g1Str);

    Fr fg, fh;
    fg.setByCSPRNG();
    fh.setByCSPRNG();
    auto g = P * fg;
    auto h = P * fh;

    cout << g << endl;
    cout << g * (-1) << endl;
    cout << g + g * (-1) << endl;
    cout << h << endl;

    uint deg = stoi(argv[1]);
    auto pp = trustedSetup(g, h, deg);
    
    Polynomial<Fr> p1 = randomPolynomial(deg);
    Polynomial<Fr> p2 = randomPolynomial(deg);

    auto com = commitPoly(p1, p2, pp.first, pp.second);
    cout << com << endl;

    Fr y;
    Fr x = 114514;
    provePolyEval(p1, x, y, pp.first);

    return 0;
}