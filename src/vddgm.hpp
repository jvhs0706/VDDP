#ifndef VDDGM_HPP
#define VDDGM_HPP

#include "vddlm.hpp"

vector<vector<bool>> BitDecomp(const vector <uint>& u, const LegendrePRNGPubParam& pp, 
    Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm,
    const Polynomial& F_in, const Polynomial& R_in, const G1& com_in,
    vector<Polynomial>& F_out, vector<Polynomial>& R_out, vector<G1>& com_out);

vector<vector<bool>> BitDecomp(const vector <int>& u, const LegendrePRNGPubParam& pp, 
    Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm,
    const Polynomial& F_in, const Polynomial& R_in, const G1& com_in,
    vector<Polynomial>& F_out, vector<Polynomial>& R_out, vector<G1>& com_out);

vector<uint> Abs(const vector<int>& s, const LegendrePRNGPubParam& pp,
    Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm,
    const Polynomial& F_in, const Polynomial& R_in, const G1& com_in,
    Polynomial& F_out, Polynomial& R_out, G1& com_out);

vector<bool> Bexp(vector<uint>u, double r, uint prob_nbit,
    const LegendrePRNGPubParam& pp, 
    Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm, 
    const Polynomial& F_in, const Polynomial& R_in, const G1& com_in,
    Polynomial& F_out, Polynomial& R_out, G1& com_out);

vector<int> DLap4DGauss(double t, uint prob_nbit,
    const LegendrePRNGPubParam& pp, 
    Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm, 
    Polynomial& F_out, Polynomial& R_out, G1& com_out);

pair<vector<int>,vector<bool>> DGaussIter(double sigma, uint prob_nbit,
    const LegendrePRNGPubParam& pp, 
    Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm, 
    Polynomial& F_out, Polynomial& R_out, G1& com_out,
    Polynomial& F_acc, Polynomial& R_acc, G1& com_acc);

vector<int> DGauss(double sigma, uint prob_nbit,
    const LegendrePRNGPubParam& pp, 
    Timer& comp_timer, Timer& ptimer, Timer& vtimer, uint& comm, 
    Polynomial& F_out, Polynomial& R_out, G1& com_out);

#endif // VDDGM_HPP