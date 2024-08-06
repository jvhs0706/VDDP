#include "utils.hpp"
using namespace std;
using namespace mcl::bn;

Polynomial randomPolynomial(uint deg) {
    vector<Fr> coefs;
    for (int i = 0; i < deg + 1; i++) {
        Fr c;
        c.setByCSPRNG();
        coefs.push_back(c);
    }
    return coefs;
}

Fr getRootOfUnity(uint logN)
{
    if (logN == 0) return 1; // X - 1
    else {
        Fr omega = -1; // X^2 - 1
        for (uint i = 2; i <= logN; ++i)
        {
            auto b = Fr::squareRoot(omega, omega);
            assert(b);
        }
        return omega;
    }
}

uint ceilLog2(uint n)
{
    uint logN = 0;
    while ((1U << logN) < n) logN++;
    return logN;
}

vector<Fr> ntt_helper(const vector<Fr>& n, const Fr& omega)
{
    assert(n.size() == 1U << ceilLog2(n.size())); // n must be a power of 2
    if (n.size() == 1) return n;

    vector<Fr> even, odd;
    for (uint i = 0; i < n.size(); i += 2)
    {
        even.push_back(n[i]);
        odd.push_back(n[i + 1]);
    }

    auto evenTransformed = ntt_helper(even, omega * omega);
    auto oddTransformed = ntt_helper(odd, omega * omega);

    vector<Fr> transformed(n.size());
    Fr omegaPower = 1;
    for (uint i = 0; i < n.size() / 2; i++)
    {
        transformed[i] = evenTransformed[i] + omegaPower * oddTransformed[i];
        transformed[i + n.size() / 2] = evenTransformed[i] - omegaPower * oddTransformed[i];
        omegaPower /= omega;
    }
    
    return transformed;
}

Polynomial ntt(const vector<Fr>& n, Fr& omega, bool inverse)
{
    omega = getRootOfUnity(ceilLog2(n.size()));
    auto transformed = ntt_helper(n, omega);
    Polynomial out(transformed);
    if (inverse) out /= Fr(n.size());
    return out;
}

Polynomial ntt_given_omega(const vector<Fr>& a, const Fr& omega, bool inverse)
{
    auto transformed = ntt_helper(a, omega);
    Polynomial out(transformed);
    if (inverse) out /= Fr(a.size());
    return out;
}