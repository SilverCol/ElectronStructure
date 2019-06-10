#include <iostream>
#include <fstream>
#include <gsl/gsl_sf_laguerre.h>
#include "Helium.h"

static const double R = 1000;
static const size_t N = 1000000;
static const double h = R/N;

void writeBinary(std::vector<double>& data, const std::string& file)
{
    std::ofstream output(file, std::ios::binary);
    for (double& x : data)
    {
        output.write(reinterpret_cast<char*>(&x), sizeof(x));
    }
    output.close();
}

int main()
{
    std::vector<double> domain(N);
    for (size_t r = 0; r < N; ++r) domain[r] = h*(r + 1);

    std::vector<double> init(N);
    for (size_t r = 0; r < init.size(); ++r)
    {
        init[r] = domain[r] * std::exp(-domain[r]);
    }
    init /= norm(init, h);

    Helium atom(init, R);

    // std::vector<double> pot = atom.electrostatic();
    // writeBinary(pot, "../data/potential.bin");
    std::vector<double> psi = atom.psi();
    writeBinary(psi, "../data/state.bin");

    std::cout << "Results: " << std::endl;
    std::cout << atom.epsilon() << std::endl;
    std::cout << atom.energy() << std::endl;
    return 0;
}