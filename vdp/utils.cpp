#include "utils.hpp"
#include <fstream>
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
            assert(Fr::squareRoot(omega, omega));
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

Polynomial ntt_vec_to_poly(const vector<Fr>& n, Fr& omega)
{
    omega = getRootOfUnity(ceilLog2(n.size()));
    auto transformed = ntt_helper(n, omega);
    Polynomial out(transformed);
    out /= Fr(n.size());
    return out;
}

Polynomial ntt_vec_to_poly_given_omega(const vector<Fr>& a, const Fr& omega)
{
    auto transformed = ntt_helper(a, omega);
    Polynomial out(transformed);
    out /= Fr(a.size());
    return out;
}

// Below are for file I/O

void savebin(const string& filename, const void* data, uint size)
{
    ofstream out(filename, ios::binary);
    out.write((char*)data, size);
    out.close();
}

uint findsize(const string& filename)
{
    ifstream in(filename, ios::binary | ios::ate);
    return in.tellg();
}

void loadbin(const string& filename, void* data, uint size)
{
    ifstream in(filename, ios::binary);
    in.read((char*)data, size);
    in.close();
}