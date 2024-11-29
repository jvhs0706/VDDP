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

    uint dim = stoi(argv[1]);
    string i_file = argv[2];
    string o_file = argv[3];
    string config_file = argv[4];

    Timer setup_timer, computing_timer;
    Timer ptimer, vtimer;
    uint comm = 0;

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

    computing_timer.start();
    Polynomial F_in = ntt_vec_to_poly_given_omega(vector<Fr>(counts.begin(), counts.end()), pp.omega_gen);
    Polynomial R_in = randomPolynomial(F_in.getDegree());
    G1 com_in = commitPoly(F_in, R_in, pp.pp.gVec, pp.pp.hVec);
    computing_timer.stop();

    Polynomial F_noise, R_noise;
    G1 com_noise;
    vector<int> noise = DiscreteLaplacianNew(z_config, g_config, pp, computing_timer, ptimer, vtimer, comm, F_noise, R_noise, com_noise);

    ptimer.start();
    Polynomial F_out = F_in + F_noise;
    Polynomial R_out = R_in + R_noise;
    vtimer.start();
    G1 com_out = com_in + com_noise;
    vtimer.stop();
    vector<int> res(dim);
    for (uint i = 0; i < dim; ++i) {
        res[i] = counts[i] + noise[i];
    }
    ptimer.stop();

    comm += (sizeof(Fr) * dim * 2);

    // save the output, binary file
    savebin(o_file, res.data(), dim * sizeof(int));

    // check output matches
    vtimer.start();
    assert(F_out == ntt_vec_to_poly_given_omega(vector<Fr>(res.begin(), res.end()), pp.omega_gen));
    vtimer.stop();

    cout << "Setting time: " << setup_timer.getTotalTime() << " s\n";
    cout << "Computing time: " << computing_timer.getTotalTime() << " s\n";
    cout << "Proving time: " << ptimer.getTotalTime() << " s\n";
    cout << "Verifying time: " << vtimer.getTotalTime() << " s\n";
    cout << "Communication: " << comm << " B\n";

    return 0;
}