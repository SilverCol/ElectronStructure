#include <iostream>
#include <fstream>
#include "Helium.h"

static const double R = 10;
static const size_t N = 100000;
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
    std::vector<double> init(N);
    for (size_t r = 0; r < init.size(); ++r) init[r] = (h*(r+1)) * std::exp(-h*(r+1));
    init /= norm(init, h);

    Helium atom(init, R);

    std::vector<double> data = atom.electrostatic();
    writeBinary(data, ".o /data/potential.bin");
    std::cout << atom.energy() << std::endl;
    std::cout << atom.norm() << std::endl;
    return 0;
}