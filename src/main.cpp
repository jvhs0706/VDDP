#include <mcl/bls12_381.hpp>
#include <iostream>
using namespace std;
using namespace mcl::bn;

#include "polynomial.hpp"
#include "timer.hpp"

int main() {
    initPairing(mcl::BLS12_381);

    Polynomial<Fr> p ({1,2,3,4,5});

    cout << p(1) << endl;
    cout << p(10000) << endl;

    return 0;
}