//
// Created by mitja on 10.6.2019.
//

#include <iostream>
#include "Helium.h"

namespace
{
    const double tolerance = 1e-4;
    const double C1 = .0;
    const double C2 = 11.4 * std::pow(8 * 3.141592653589793 / 3, 1.0/3.0);
    const double C = std::pow(3 / (2 * std::pow(3.141592653589793, 2)), 1.0/3.0);
}

Helium::Helium(const std::vector<double>& u, double R):
m_u(u),
m_r(u.size()),
m_V(u.size()),
m_f(u.size()),
m_h(R/u.size()),
m_e(0.0),
m_N(u.size())
{
    constructSpace();
}

std::vector<double> Helium::electrostatic()
{
    if (std::abs(norm() - 1.0) > tolerance)
    {
        throw std::runtime_error("State not normalized.");
    }

    for (size_t n = 0; n < m_N; ++n)
    {
        m_f[n] = - std::pow(m_h*m_u[n], 2) / (12 * m_r[n]);
    }

    std::vector<double> U(m_N, 0.0);
    U[1] = m_h + (7*m_f[0] + 6*m_f[1] - m_f[2]) / 2;

    for (size_t n = 2; n < m_N; ++n)
    {
        U[n] = 2*U[n-1] - U[n-2] + m_f[n] + 10*m_f[n-1] + m_f[n-2];
    }

    double k = (1 - U.back()) / m_r.back();
    for (size_t r = 0; r < m_N; ++r)
    {
        U[r] += k * m_r[r];
    }

    return U;
}

std::vector<double> Helium::correlatic()
{
    std::vector<double> v(m_N);
    std::transform(m_u.begin(), m_u.end(), m_r.begin(), v.begin(), [](double u, double r)
    {
        return - C * std::pow(std::abs(u / r), 2.0/3.0)
                + C1 * std::log(1 + C2 * std::pow(std::abs(u / r), 2.0/3.0));
    });

    return v;
}

void Helium::LDA_DFT()
{
    std::cout << "Running LDA DFT." << std::endl;
    uint32_t count = 0;
    while(true)
    {
        ++ count;
        std::cout << "Iteration no. " << count <<  " | Now at " << m_e << '\r' << std::flush;

        //std::vector<double> comparison(m_u);
        double comparison = m_e;
        std::vector<double> U = electrostatic();
        std::vector<double> V = correlatic();

        for (size_t n = 0; n < m_N; ++n)
        {
            m_V[n] = (2*U[n] - 2) / m_r[n] + V[n];
        }

        updateDensity();
        if (std::abs(m_e - comparison) < tolerance)
        {
            std::cout << "Finished in " << count << " iterations." << std::endl;
            return;
        }
    }
}

double Helium::energy()
{
    std::vector<double> U = electrostatic();
    std::vector<double> V = correlatic();

    double energy = 0.0;
    for (size_t n = 0; n < m_N; ++n)
    {
        energy -= 2 * U[n] * std::pow(m_u[n], 2) / m_r[n];
        energy -= V[n] * std::pow(m_u[n], 2) / 2;
    }
    energy *= m_h;
    return 2*m_e + energy;
}

std::vector<double> Helium::psi()
{
    std::vector<double> psi(m_N);
    std::transform(m_u.begin(), m_u.end(), m_r.begin(), psi.begin(), std::multiplies<>());
    return psi;
}
