#include <mcl/bls12_381.hpp>
#include <iostream>
using namespace std;
using namespace mcl::bn;

#include "utils.hpp"
#include "polynomial.hpp"
#include "proof.hpp"
#include "timer.hpp"
#include "vrr.hpp"

int main(int argc, char** argv) {
    initPairing(mcl::BLS12_381);
    // initialize the random seed
    srand(time(NULL));

    uint num_class = 8;
    uint omega_size = 64;
    uint B = 3;
    uint A = omega_size - (num_class - 1) * B;
    auto vrrpp = VRRTrustedSetup(num_class, A, B);

    // random integer in [0, num_class)
    uint ix = rand() % num_class;
    // random inter in [0, omega_size)
    uint is = rand() % omega_size;

    Fr rx, rs;
    rx.setByCSPRNG();
    rs.setByCSPRNG();

    auto com = VRRCommit(ix, is, rx, rs, vrrpp);
    auto& comx = com.first;
    auto& coms = com.second;

    uint ir = rand() % omega_size;
    auto y = VRRCompute(ix, is, ir, vrrpp);

    cout << y << endl;

    Timer ptimer, vtimer;

    bool res = VRR(ix, is, rx, rs, comx, coms, y, ir, vrrpp, ptimer, vtimer);


    return 0;
}