#include "prg.hpp"

int LegendreSymbol(const Fr& x)
{
    if (x.isZero()) return 0;
    else 
    {
        Fr temp;
        bool b = Fr::squareRoot(temp, x);
        if (b) return 1;
        else return -1;
    }
}

bool LegendrePRF(const Fr& x, const Fr& key)
{
    return LegendreSymbol(x + key) >= 0;
}