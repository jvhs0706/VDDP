#ifndef VDDLM_HPP
#define VDDLM_HPP

#include <iostream>
#include <vector>
#include "proof.hpp"
#include "prng.hpp"

using namespace std;

vector<bool> probToBits(double p, uint n);

void readVDDLMConfig(const string& config_file, vector<bool>& z_config, vector<vector<bool>>& g_config);

void readVDDLMConfig(const string& config_file, vector<bool>& z_config, vector<vector<bool>>& g_config);

vector<bool> Bernoulli(vector<bool> p, const LegendrePRNGPubParam& pp, Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm, Polynomial& F_out, Polynomial& R_out, G1& com_out);

vector<uint> Geometric(double p, uint log_range, uint prec, const LegendrePRNGPubParam& pp, Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm, Polynomial& F_out, Polynomial& R_out, G1& com_out);

vector<int> DiscreteLaplacian(double t, uint log_range, uint prec, const LegendrePRNGPubParam& pp, Timer& comp_timer, Timer& ptimer, Timer& vtimer, Polynomial& F_out, Polynomial& R_out, G1& com_out);

vector<uint> GeometricNew(const vector<vector<bool>>& config, const LegendrePRNGPubParam& pp, Timer& comp_timer, Timer& ptimer, Timer& vtimer, Polynomial& F_out, Polynomial& R_out, G1& com_out);

vector<int> DiscreteLaplacianNew(const vector<bool>& z_config, const vector<vector<bool>>& g_config, const LegendrePRNGPubParam& pp, Timer& comp_timer, Timer& ptimer, Timer& vtimer, Polynomial& F_out, Polynomial& R_out, G1& com_out);

#endif // VDDLM_HPP