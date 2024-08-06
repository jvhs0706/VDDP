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

    // we need to find a non-square first
    Fr temp = 0;
    auto b = Fr::squareRoot(temp, temp);
    cout << b << " " << temp << endl;

    Fr key;
    key.setByCSPRNG();
    vector<Fr> rt_vec;
    uint len = stoi(argv[1]);
    vector<bool> res = LegendrePRNG(key, len, rt_vec);

    for (uint i = 0; i < len; ++ i)
    {
        cout << i + key << " " << rt_vec[i] * rt_vec[i] << " " << rt_vec[i] * rt_vec[i] / LEGENDRE_PRNG_NON_SQ << endl;
    }


    // uint len = stoi(argv[1]);
    // double t = stod(argv[2]);
    // uint log_range = stoi(argv[3]);
    // uint prec = stoi(argv[4]);

    // // test the LegendrePRF function

    // vector<vector<bool>> r1;
    // vector<vector<vector<bool>>> r2(log_range);

    // Timer timer;

    // timer.start();
    // Fr key;

    // for (uint j = 0; j < prec; ++ j)
    // {
    //     key.setByCSPRNG();
    //     r1.push_back(LegendrePRNG(key, len));
    // }

    // for (uint i = 0; i < log_range; ++ i)
    // {
    //     for (uint j = 0; j < prec; ++ j)
    //     {
    //         key.setByCSPRNG();
    //         r2[i].push_back(LegendrePRNG(key, len));
    //     }
    // }

    // key.setByCSPRNG();
    // vector<bool> r3 = LegendrePRNG(key, len);

    // timer.stop();
    // cout << "PRNG time = " << timer.getTotalTime() << " s" << endl;
    // timer.reset();

    // timer.start();
    // auto res = sampleLaplacian(r1, r2, r3, len, t, log_range, prec);
    // timer.stop();
    // cout << "Laplacian time = " << timer.getTotalTime() << " s" << endl;

    // for (uint i = 0; i < len; ++ i)
    // {
    //     cout << res[i] << " ";
    // }

    return 0;
}