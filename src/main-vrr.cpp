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

    uint num_class = stoi(argv[1]);
    uint omega_size = stoi(argv[2]);
    uint A = stoi(argv[3]);
    uint B = stoi(argv[4]);

    assert (A + (num_class - 1)*B == omega_size);
    
    Timer setup_timer, committing_timer, computing_timer;

    setup_timer.start();
    auto vrrpp = VRRTrustedSetup(num_class, A, B);
    setup_timer.stop();

    // random integer in [0, num_class)
    uint ix = rand() % num_class;
    // random inter in [0, omega_size)
    uint is = rand() % omega_size;

    Fr rx, rs;
    rx.setByCSPRNG();
    rs.setByCSPRNG();

    committing_timer.start();
    auto com = VRRCommit(ix, is, rx, rs, vrrpp);
    committing_timer.stop();

    auto& comx = com.first;
    auto& coms = com.second;

    computing_timer.start();
    uint ir = rand() % omega_size;
    auto y = VRRCompute(ix, is, ir, vrrpp);
    computing_timer.stop();

    Timer ptimer, vtimer;

    uint comm;

    assert(VRR(ix, is, rx, rs, comx, coms, y, ir, vrrpp, ptimer, vtimer, comm));

    cout << "Setting time: " << setup_timer.getTotalTime() << " s\n";
    cout << "Committing time: " << committing_timer.getTotalTime() << " s\n";
    cout << "Computing time: " << computing_timer.getTotalTime() << " s\n";
    cout << "Proving time: " << ptimer.getTotalTime() << " s\n";
    cout << "Verifying time: " << vtimer.getTotalTime() << " s\n";
    cout << "Communication: " << comm << " bytes\n";


    return 0;
}