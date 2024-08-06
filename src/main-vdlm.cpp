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

    Polynomial F ({1, 9, 2, 6, 0, 8, 1, 7});
    VanishingPolynomial G (4);
    Polynomial H ({1, 1, 4, 5, 1, 4});
    Polynomial R = F * H + G;

    cout << R / F << endl;
    cout << R % F << endl;

    return 0;
}