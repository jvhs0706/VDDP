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
    double prob = stod(argv[2]);
    uint prec = stoi(argv[3]);

    // test the LegendrePRF function

    vector<vector<bool>> r_vec;

    Timer timer;

    timer.start();
    Fr key;
    for (uint i = 0; i < prec; ++ i)
    {
        key.setByCSPRNG();
        r_vec.push_back(LegendrePRNG(key, len));
    }
    timer.stop();
    cout << "PRNG time = " << timer.getTotalTime() << " s" << endl;
    timer.reset();

    
    timer.start();
    auto res = sampleBernoulli(r_vec, len, prob, prec);
    timer.stop();
    cout << "bernoulli time = " << timer.getTotalTime() << " s" << endl;

    int res_sum = 0;
    for (uint i = 0; i < len; ++ i)
    {
        res_sum += res[i];
    }


    return 0;
}