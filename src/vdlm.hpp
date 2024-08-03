#ifndef VDLP_HPP
#define VDLP_HPP

#include <iostream>
#include <vector>

using namespace std;

vector<bool> probToBits(double p, uint n);

// if a[i] is true, res[i] = b[i], otherwise res[i] = c[i]
// in arithmetic: res[i] = a[i] * b[i] + (1 - a[i]) * c[i] = a[i] * (b[i] - c[i]) + c[i]
// vector<bool> mux(const vector<bool>& a, const vector<bool>& b, const vector<bool>& c);

vector<bool> sampleBernoulli(const vector<vector<bool>>& r_vec, uint len, double p, uint prec);

vector<uint> sampleGeometric(const vector<vector<vector<bool>>>& r_vec, uint len, double p, uint log_range, uint prec);

vector<int> sampleLaplacian(const vector<vector<bool>>& r1, const vector<vector<vector<bool>>>& r2, const vector<bool>& r3,
    uint len, double t, uint log_range, uint prec);

#endif // VDLP_HPP