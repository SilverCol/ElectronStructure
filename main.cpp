#include <iostream>
#include "Helium.h"

static const double R = 10;
static const size_t N = 100000;
static const double h = R/N;

int main()
{
    std::vector<double> init(N);
    for (size_t r = 0; r < init.size(); ++r) init[r] = (h*(r+1)) * std::exp(-h*(r+1));
    init /= norm(init, h);

    Helium atom(init, R);

    std::cout << atom.energy() << std::endl;
    std::cout << atom.norm() << std::endl;
    return 0;
}