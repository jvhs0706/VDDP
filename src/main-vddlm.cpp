#include <mcl/bls12_381.hpp>
#include <iostream>
using namespace std;
using namespace mcl::bn;

#include "utils.hpp"
#include "polynomial.hpp"
#include "proof.hpp"
#include "timer.hpp"
#include "prng.hpp"
#include "vddlm.hpp"

int main(int argc, char** argv) {
    initPairing(mcl::BLS12_381);
    // initialize the random seed
    srand(time(NULL));

    uint len = stoi(argv[1]);
    auto pp = LegendrePRNGTrustedSetup(len);

    float p = stod(argv[2]);
    uint prec = stoi(argv[3]);
    uint log_range = stoi(argv[4]);

    Timer setup_timer, computing_timer;
    Timer ptimer, vtimer;

    Polynomial F_out, R_out;
    G1 com_out;

    Polynomial F, R;
    G1 com;

    auto p_bin = probToBits(p, prec);
    // print out p_bin
    for (auto b : p_bin) cout << b;
    cout << endl;

    auto g = Geometric(p, log_range, prec, pp, computing_timer, ptimer, vtimer, F, R, com_out);
    uint sum = 0;
    for (const auto & i : g) {
        sum += i;
    }
    // print the mean of g
    cout << "mean = " << (double)sum / len << endl;
    cout << computing_timer.getTotalTime() << " " << ptimer.getTotalTime() << " " << vtimer.getTotalTime() << endl;
    return 0;

    

    return 0;
}