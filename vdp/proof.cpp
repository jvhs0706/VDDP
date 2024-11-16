#include "proof.hpp"
#include "utils.hpp"
using namespace std;
using namespace mcl::bn;

PubParam trustedSetup(uint maxDeg)
{
    G1 P;
    const char *g1Str = "1 0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb 0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1";
    P.setStr(g1Str);

    G2 Q;
    const char *g2Str = "1 0x24aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8 0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e 0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801 0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be";
    Q.setStr(g2Str);

    Fr fg, fh, fg2;
    fg.setByCSPRNG();
    fh.setByCSPRNG();
    fg2.setByCSPRNG();

    auto g = P * fg;
    auto h = P * fh;
    auto g2 = Q * fg2;

    Fr tau;
    tau.setByCSPRNG();
    Fr temp = 1;
    G2 g2_tau = g2 * tau;
    vector<G1> gVec(maxDeg + 1), hVec(maxDeg + 1);
    for (uint i = 0; i <= maxDeg; i++)
    {
        gVec[i] = g * temp;
        hVec[i] = h * temp;
        temp *= tau;
    }
    return {gVec, hVec, g2, g2_tau};
}

G1 commitPoly(const Polynomial& F, const std::vector<G1>& gVec)
{
    if (F.getDegree() == -1)
    {
        return gVec[0] * 0;
    }
    assert(F.getDegree() < gVec.size());

    G1 result = gVec[0] * F[0];
    for (uint i = 1; i <= F.getDegree(); i++)
    {
        result += gVec[i] * F[i];
    }
    return result;

}

G1 commitPoly(const Polynomial& F, const Polynomial& R, const std::vector<G1>& gVec, const std::vector<G1>& hVec)
{
    if (R.getDegree() == -1)
    {
        return commitPoly(F, gVec);
    }
    if (F.getDegree() == -1)
    {
        return commitPoly(R, hVec);
    }
    assert(F.getDegree() < gVec.size());
    assert(R.getDegree() >= 0);
    assert(R.getDegree() < hVec.size());

    G1 result = gVec[0] * F[0] + hVec[0] * R[0];
    for (uint i = 1; i <= F.getDegree(); i++)
    {
        result += gVec[i] * F[i];
    }
    for (uint i = 1; i <= R.getDegree(); i++)
    {
        result += hVec[i] * R[i];
    }
    return result;
}


mcl::bn::G1 provePolyEval(const Polynomial& F, 
    const mcl::bn::Fr& x, mcl::bn::Fr& y, const std::vector<mcl::bn::G1>& gVec)
{
    y = F(x);
    auto diff = Polynomial({-x, 1});
    return commitPoly((F - y) / diff, gVec);
}

mcl::bn::G1 provePolyEval(const Polynomial& F, const Polynomial& R, 
    const mcl::bn::Fr& x, 
    mcl::bn::Fr& yF, mcl::bn::Fr& yR,
    const std::vector<mcl::bn::G1>& gVec, const std::vector<mcl::bn::G1>& hVec)
{
    yF = F(x);
    yR = R(x);
    auto diff = Polynomial({-x, 1});
    return commitPoly((F - yF) / diff, (R - yR) / diff, gVec, hVec);
}

bool verifyPolyEval(const Fr& yF, const Fr& yR, const G1& proof, const G1& commitment, const Fr& x, const G1& g, const G1& h, const G2& g2, const G2& g2_tau)
{
    auto g2_diff = g2_tau - g2 * x;

    GT l, r;
    pairing(l, commitment - g * yF - h * yR, g2);
    pairing(r, proof, g2_diff);

    return l == r;
}

bool verifyPolyEval(const Fr& y, const G1& proof, const G1& commitment, const Fr& x, const G1& g, const G2& g2, const G2& g2_tau)
{
    auto g2_diff = g2_tau - g2 * x;

    GT l, r;
    pairing(l, commitment - g * y, g2);
    pairing(r, proof, g2_diff);

    return l == r;
}

bool SecretEval(const Fr& y, const Fr& ry, const Fr& x, const Fr& rx, const Polynomial& F, 
    const G1& com_y, const G1& com_x, const G1& com_F, 
    const PubParam& pp, Timer& ptimer, Timer& vtimer)
{
    bool accepted = true;

    const auto& gVec = pp.gVec;
    const auto& hVec = pp.hVec;
    const auto& g2 = pp.g2;
    const auto& g2_tau = pp.g2_tau;
    const auto& g = gVec[0];
    const auto& h = hVec[0];

    assert(y==F(x));
    assert(com_x == g * x + h * rx);
    assert(com_y == g * y + h * ry);

    ptimer.start();
    auto diff = Polynomial({x, -1});
    auto F_ = (-F + y) / diff;
    auto R_F_ = randomPolynomial(gVec.size() - 1);
    auto com_F_ = commitPoly(F_, R_F_, gVec, hVec);
    ptimer.stop();

    vtimer.start();
    Fr u;
    u.setByCSPRNG();
    vtimer.stop();

    ptimer.start();
    auto z = F(u), z_ = F_(u);
    Fr r_z_;
    r_z_.setByCSPRNG();
    auto com_z_ = g * z_ + h * r_z_;
    ptimer.stop();

    accepted &= Prod(y - z, ry, x - u, rx, z_, r_z_, com_y - g * z, com_x - g * u, com_z_, g, h, ptimer, vtimer);

    ptimer.start();
    auto proof = provePolyEval(F, u, z, gVec);
    ptimer.stop();

    vtimer.start();
    accepted &= verifyPolyEval(z, proof, com_F, u, g, g2, g2_tau);
    vtimer.stop();

    ptimer.start();
    auto r_ = R_F_(u) - r_z_;
    Fr yF_, yR_;
    auto proof_ = provePolyEval(F_ - z_, R_F_ - r_z_, u, yF_, yR_, gVec, hVec);
    ptimer.stop();

    vtimer.start();
    accepted &= (yF_ == 0);
    accepted &= verifyPolyEval(yF_, yR_, proof_, com_F_ - com_z_, u, g, h, g2, g2_tau);
    vtimer.stop();

    return accepted;
}

bool EvalSecret(const Fr& y, const Fr& ry, const Fr& x, // secret y, public x
    const Polynomial& F, const Polynomial& R, // secret F, R
    const G1& com_y, const G1& com_F, 
    const PubParam& pp, Timer& ptimer, Timer& vtimer)
{
    bool accepted = true;

    const auto& gVec = pp.gVec;
    const auto& hVec = pp.hVec;
    const auto& g2 = pp.g2;
    const auto& g2_tau = pp.g2_tau;
    const auto& g = gVec[0];
    const auto& h = hVec[0];

    assert(y==F(x));
    assert(com_y == g * y + h * ry);

    ptimer.start();
    auto F_ = F - y;
    auto R_ = R - ry;
    vtimer.start();
    auto com_F_ = com_F - com_y;
    vtimer.stop();
    ptimer.stop();

    Fr y_, r_;
    ptimer.start();
    auto proof = provePolyEval(F_, R_, x, y_, r_, gVec, hVec);
    ptimer.stop();

    vtimer.start();
    accepted &= y_.isZero();
    accepted &= verifyPolyEval(y_, r_, proof, com_F_, x, g, h, g2, g2_tau);
    vtimer.stop();

    return accepted;
}

bool Prod(const Fr& z, const Fr& r_z, const Fr& x, const Fr& r_x, const Fr& y, const Fr& r_y,
    const G1& com_z, const G1& com_x, const G1& com_y, const G1& g, const G1& h,
    Timer& ptimer, Timer& vtimer)
    
{
    bool accepted = true;
    ptimer.start();
    Fr b1, b2, b3, b4, b5;
    b1.setByCSPRNG();
    b2.setByCSPRNG();
    b3.setByCSPRNG();
    b4.setByCSPRNG();
    b5.setByCSPRNG();
    auto alpha = g * b1 + h * b2;
    auto beta = g * b3 + h * b4;
    auto delta = com_x * b3 + h * b5;
    ptimer.stop();

    vtimer.start();
    Fr c;
    c.setByCSPRNG();
    vtimer.stop();

    ptimer.start();
    auto z1 = b1 + c * x;
    auto z2 = b2 + c * r_x;
    auto z3 = b3 + c * y;
    auto z4 = b4 + c * r_y;
    auto z5 = b5 + c * (r_z - r_x * y);
    ptimer.stop();

    vtimer.start();
    accepted &= (alpha + com_x * c == g * z1 + h * z2);
    accepted &= (beta + com_y * c == g * z3 + h * z4);
    accepted &= (delta + com_z * c == com_x * z3 + h * z5);
    vtimer.stop();

    return accepted;
}

bool Equal(const Fr& x, const Fr& r_x, const Fr& y, const Fr& r_y, 
    const G1& com_x, const G1& com_y, const G1& g, const G1& h,
    Timer& ptimer, Timer& vtimer)
{
    bool accepted = true;

    ptimer.start();
    Fr s;
    s.setByCSPRNG();
    G1 a = h * s;
    ptimer.stop();

    vtimer.start();
    Fr c;
    c.setByCSPRNG();
    vtimer.stop();

    ptimer.start();
    auto z = c * (r_x - r_y) + s;
    ptimer.stop();

    vtimer.start();
    accepted &= (h * z == (com_x - com_y) * c + a);
    vtimer.stop();

    return accepted;
}

bool Binary(const Polynomial& F, const Polynomial& R, 
    const uint len, const Fr& omega, 
    const G1& com_F, const PubParam& pp, 
    Timer& ptimer, Timer& vtimer)
{
    bool accepted = true;
    ptimer.start();
    Polynomial zero_poly = F * F - F;
    VanishingPolynomial V(len);
    Polynomial F_quot, F_rem;
    zero_poly.divide(V, F_quot, F_rem);
    assert (F_rem.getDegree() == -1); // It's the zero polynomial
    auto R_quot = randomPolynomial(len);
    auto com_F_quot = commitPoly(F_quot, R_quot, pp.gVec, pp.hVec);
    ptimer.stop();

    vtimer.start();
    Fr u;
    u.setByCSPRNG();
    vtimer.stop();

    ptimer.start();
    const auto x = F(u), z_ = F_quot(u), y = V(u); // z = x * y
    Fr r_x, r_z_;
    r_x.setByCSPRNG();
    r_z_.setByCSPRNG();
    
    const G1& g = pp.gVec[0];
    const G1& h = pp.hVec[0];
    const auto com_x = g * x + h * r_x;
    const auto com_z_ = g * z_ + h * r_z_;
    const auto z = z_ * y + x;
    const auto r_z = r_z_ * y + r_x;
    const auto com_z = com_z_ * y + com_x;
    
    ptimer.stop();

    accepted &= EvalSecret(x, r_x, u, F, R, com_x, com_F, pp, ptimer, vtimer);
    accepted &= EvalSecret(z_, r_z_, u, F_quot, R_quot, com_z_, com_F_quot, pp, ptimer, vtimer);
    accepted &= Prod(z, r_z, x, r_x, x, r_x, com_z, com_x, com_x, g, h, ptimer, vtimer);

    return accepted;
}

// c = a * b
bool Hadamard(const Polynomial& Fc, const Polynomial& Rc, 
    const Polynomial& Fa, const Polynomial& Ra,
    const Polynomial& Fb, const Polynomial& Rb,
    const uint len, const Fr& omega,
    const G1& com_Fc, const G1& com_Fa, const G1& com_Fb,
    const PubParam& pp,
    Timer& ptimer, Timer& vtimer)    
{
    bool accepted = true;
    ptimer.start();
    Polynomial zero_poly = Fa * Fb - Fc;
    VanishingPolynomial V(len);
    Polynomial F_quot, F_rem;
    zero_poly.divide(V, F_quot, F_rem);
    assert (F_rem.getDegree() == -1); // It's the zero polynomial
    auto R_quot = randomPolynomial(len);
    auto com_F_quot = commitPoly(F_quot, R_quot, pp.gVec, pp.hVec);
    ptimer.stop();

    vtimer.start();
    Fr u;
    u.setByCSPRNG();
    const auto v = V(u);
    vtimer.stop();

    ptimer.start();
    const auto ya = Fa(u), yb = Fb(u), yc = Fc(u), y_quot = F_quot(u); // z = x * y
    Fr ra, rb, rc, r_quot;
    ra.setByCSPRNG();
    rb.setByCSPRNG();
    rc.setByCSPRNG();
    r_quot.setByCSPRNG();

    const G1& g = pp.gVec[0];
    const G1& h = pp.hVec[0];
    const auto com_a = g * ya + h * ra;
    const auto com_b = g * yb + h * rb;
    const auto com_c = g * yc + h * rc;
    const auto com_quot = g * y_quot + h * r_quot;
    ptimer.stop();

    accepted &= EvalSecret(ya, ra, u, Fa, Ra, com_a, com_Fa, pp, ptimer, vtimer);
    accepted &= EvalSecret(yb, rb, u, Fb, Rb, com_b, com_Fb, pp, ptimer, vtimer);
    accepted &= EvalSecret(yc, rc, u, Fc, Rc, com_c, com_Fc, pp, ptimer, vtimer);
    accepted &= EvalSecret(y_quot, r_quot, u, F_quot, R_quot, com_quot, com_F_quot, pp, ptimer, vtimer);
    accepted &= Prod(v * y_quot + yc, v * r_quot + rc, yb, rb, ya, ra, com_quot * v + com_c, com_b, com_a, g, h, ptimer, vtimer);

    return accepted;
}

//out - b = s * (a - b)
bool Mux(const Polynomial& Fs, const Polynomial& Rs, // selector
    const Polynomial& Fa, const Polynomial& Ra,
    const Polynomial& Fb, const Polynomial& Rb,
    const Polynomial& Fout, const Polynomial& Rout,
    const uint len, const Fr& omega,
    const G1& com_Fs, const G1& com_Fa, const G1& com_Fb, const G1& com_Fout,
    const PubParam& pp,
    Timer& ptimer, Timer& vtimer)
{
    ptimer.start();
    const auto F_lhs = Fout - Fb;
    const auto R_lhs = Rout - Rb;
    const auto& F_rhs_1 = Fs;
    const auto& R_rhs_1 = Rs;
    const auto F_rhs_2 = Fa - Fb;
    const auto R_rhs_2 = Ra - Rb;
    vtimer.start();
    const auto com_F_lhs = com_Fout - com_Fb;
    const auto& com_F_rhs_1 = com_Fs;
    const auto com_F_rhs_2 = com_Fa - com_Fb;
    vtimer.stop();
    ptimer.stop();

    return Hadamard(F_lhs, R_lhs, F_rhs_1, R_rhs_1, F_rhs_2, R_rhs_2, len, omega, com_F_lhs, com_F_rhs_1, com_F_rhs_2, pp, ptimer, vtimer);
}