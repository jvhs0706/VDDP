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

    uint nser = stoi(argv[1]);
    uint dim = stoi(argv[2]);
    string i_file = argv[3];
    string o_file = argv[4];
    string config_file = argv[5];

    Timer setup_timer, computing_timer;
    Timer ptimer, vtimer;

    setup_timer.start();
    auto pp = LegendrePRNGTrustedSetup(dim);
    setup_timer.stop();

    // load the input, binary file
    uint size = findsize(i_file);
    assert(size == dim * sizeof(int));
    vector<int> counts(dim);
    loadbin(i_file, counts.data(), size);

    vector<bool> z_config;
    vector<vector<bool>> g_config;
    readVDDLMConfig(config_file, z_config, g_config);

    // create additive secret shares
    vector<vector<Fr>> shares(nser, vector<Fr>(dim));
    shares[0] = vector<Fr>(counts.begin(), counts.end());
    for (uint i = 1; i < nser; ++ i)
    {
        for (uint j = 0; j < dim; ++ j)
        {
            shares[i][j].setByCSPRNG();
            shares[0][j] -= shares[i][j];
        }
    }

    vector<int> res = counts;

    for (uint i = 0; i < nser; ++ i)
    {
        const vector<Fr>& cur_share = shares[i];
        computing_timer.start();
        Polynomial F_in = ntt_vec_to_poly_given_omega(cur_share, pp.omega_gen);
        Polynomial R_in = randomPolynomial(F_in.getDegree());
        G1 com_in = commitPoly(F_in, R_in, pp.pp.gVec, pp.pp.hVec);
        computing_timer.stop();

        Polynomial F_noise, R_noise;
        G1 com_noise;
        vector<int> noise = DiscreteLaplacianNew(z_config, g_config, pp, computing_timer, ptimer, vtimer, F_noise, R_noise, com_noise);

        ptimer.start();
        Polynomial F_out = F_in + F_noise;
        Polynomial R_out = R_in + R_noise;
        vtimer.start();
        G1 com_out = com_in + com_noise;
        vector<Fr> temp_res(dim);
        for (uint j = 0; j < dim; ++ j) {
            temp_res[j] = cur_share[j] + noise[j];
            res[j] += noise[j];
        }
        vtimer.stop();
        ptimer.stop();

        vtimer.start();
        assert(commitPoly(F_out, R_out, pp.pp.gVec, pp.pp.hVec) == com_out);
        assert(F_out == ntt_vec_to_poly_given_omega(temp_res, pp.omega_gen));
        vtimer.stop();

    }

    // save the output, binary file
    savebin(o_file, res.data(), dim * sizeof(int));

    cout << "Setting time: " << setup_timer.getTotalTime() << " s\n";
    cout << "Computing time: " << computing_timer.getTotalTime() << " s\n";
    cout << "Proving time: " << ptimer.getTotalTime() << " s\n";
    cout << "Verifying time: " << vtimer.getTotalTime() << " s\n";

    return 0;
}