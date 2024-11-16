#ifndef VDDLM_HPP
#define VDDLM_HPP

#include <iostream>
#include <vector>
#include "proof.hpp"
#include "prng.hpp"

using namespace std;

vector<bool> probToBits(double p, uint n);

vector<bool> Bernoulli(vector<bool> p, const LegendrePRNGPubParam& pp, Timer& comp_timer, Timer& ptimer, Timer& vtimer, Polynomial& F_out, Polynomial& R_out, G1& com_out);

vector<uint> Geometric(double p, uint log_range, uint prec, const LegendrePRNGPubParam& pp, Timer& comp_timer, Timer& ptimer, Timer& vtimer, Polynomial& F_out, Polynomial& R_out, G1& com_out);

vector<int> DiscreteLaplacian(double t, uint log_range, uint prec, const LegendrePRNGPubParam& pp, Timer& comp_timer, Timer& ptimer, Timer& vtimer, Polynomial& F_out, Polynomial& R_out, G1& com_out);

#endif // VDDLM_HPP