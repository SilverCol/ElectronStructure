#include <iostream>
#include "Helium.h"

static const double R = 10;
static const size_t N = 100000;

int main()
{
    std::vector<double> init(N);
    Helium atom(init, R);

    std::cout << atom.energy() << std::endl;
    return 0;
}