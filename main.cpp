#include <iostream>
#include <fstream>
#include <gsl/gsl_sf_laguerre.h>
#include "SpectralHe.h"
#include "NumerovHe.h"

static const double R = 17;
static const size_t N = 10000000;
static const double h = R/N;

void writeBinary(const std::vector<double>& data, const std::string& file)
{
    std::ofstream output(file, std::ios::binary);
    for (double x : data)
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

    NumerovHe atom(init, R);

    std::string pot = "../data/potential";
    pot.append(std::to_string((int)R));
    pot.append(".bin");
    writeBinary(atom.electrostatic(), pot);
    
    std::string state = "../data/state";
    state.append(std::to_string((int)R));
    state.append(".bin");
    writeBinary(atom.psi(), state);

    std::cout << "Results: " << std::endl;
    std::cout << atom.epsilon() << std::endl;
    std::cout << atom.energy() << std::endl;
    return 0;
}