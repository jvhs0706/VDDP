#include <mcl/bls12_381.hpp>
#include <iostream>
using namespace std;
using namespace mcl::bn;

#include "utils.hpp"
#include "polynomial.hpp"
#include "proof.hpp"
#include "timer.hpp"
#include "prg.hpp"

int main(int argc, char** argv) {
    initPairing(mcl::BLS12_381);
    // initialize the random seed
    srand(time(NULL));

    // test the LegendrePRF function

    Fr key;
    key.setByCSPRNG();

    for (uint i = 0; i < 100; ++ i)
    {
        Fr x;
        x.setByCSPRNG();
        cout << LegendrePRF(x, key) << endl;
    }

    return 0;
}