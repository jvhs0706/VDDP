#ifndef VRR_HPP
#define VRR_HPP

#include <mcl/bls12_381.hpp>
#include <vector>
#include "utils.hpp"
#include "polynomial.hpp"
#include "proof.hpp"

using namespace mcl::bn;
using namespace std;

struct VRRPubParam
{
    uint num_class;
    uint omega_size;
    PubParam pp;
    Fr x_gen;
    Fr omega_gen;
    Polynomial F;
    Polynomial F_omega;
    G1 com_F;
    G1 com_F_omega; 
};

VRRPubParam VRRTrustedSetup(uint num_class, uint A, uint B);

pair<G1, G1> VRRCommit(const Fr& ix, const uint is, const Fr& rx, const Fr& rs, const VRRPubParam& vrrpp);

Fr VRRCompute(const uint ix, const uint is, const uint ir, const VRRPubParam& vrrpp);

bool VRR(const uint ix, const uint is, const Fr& rx, const Fr& rs, 
    const G1& comx, const G1& coms, const Fr& y, const uint ir, 
    const VRRPubParam& vrrpp, Timer &ptimer, Timer &vtimer);


#endif // VRR_HPP