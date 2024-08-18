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

    float p = stod(argv[2]);
    uint prec = stoi(argv[3]);

    Timer setup_timer, committing_timer, computing_timer;
    Timer ptimer, vtimer;

    G1 com_out;

    Bernoulli(probToBits(p, prec), pp, ptimer, vtimer, com_out);

    cout << "Proving time: " << ptimer.getTotalTime() << " s\n";
    cout << "Verifying time: " << vtimer.getTotalTime() << " s\n";

    return 0;
}