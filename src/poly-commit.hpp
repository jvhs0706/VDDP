#ifndef POLY_COMMIT_HPP
#define POLY_COMMIT_HPP

#include <mcl/bls12_381.hpp>
#include <vector>
#include "polynomial.hpp"


std::pair<std::vector<mcl::bn::G1>, std::vector<mcl::bn::G1>> trustedSetup(const mcl::bn::G1& g, const mcl::bn::G1& h, const mcl::bn::G2& g2, const mcl::bn::G2& g2_tau, uint maxDeg);

mcl::bn::G1 commitPoly(const Polynomial<mcl::bn::Fr>& F, const std::vector<mcl::bn::G1>& gVec);

mcl::bn::G1 commitPoly(const Polynomial<mcl::bn::Fr>& F, const Polynomial<mcl::bn::Fr>& R, const std::vector<mcl::bn::G1>& gVec, const std::vector<mcl::bn::G1>& hVec);

mcl::bn::G1 provePolyEval(const Polynomial<mcl::bn::Fr>& F, 
    const mcl::bn::Fr& x, mcl::bn::Fr& y, const std::vector<mcl::bn::G1>& gVec);

mcl::bn::G1 provePolyEval(const Polynomial<mcl::bn::Fr>& F, const Polynomial<mcl::bn::Fr>& R, 
    const mcl::bn::Fr& x, 
    mcl::bn::Fr& yF, mcl::bn::Fr& yR,
    const std::vector<mcl::bn::G1>& gVec, const std::vector<mcl::bn::G1>& hVec);


#endif // POLY_COMMIT_HPP