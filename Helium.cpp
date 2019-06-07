//
// Created by mitja on 4.6.2019.
//

#include <iostream>
#include "Helium.h"


namespace
{
    const double epsilon = 1e-2;
    const double F = std::pow(3.141592653589793 / 6, 1.0/3.0);
    const double A = F * .0545 * std::pow(8.0 * 3.141592653589793 / 3.0, -1.0/3.0);
    const double B = 11.4 * std::pow(8.0 * 3.141592653589793 / 3.0, 1.0/3.0);
}

Helium::Helium(const std::vector<double>& u, double R) :
m_u(u),
m_V(u.size()),
m_f(u.size()),
m_h(R/u.size()),
m_E(0.0)
{
    for (size_t r = 0; r < m_V.size(); ++r) m_V[r] = - 1/((r + 1)*m_h);
    updateDensity();
    // LDA_DFT();
}

void Helium::shoot(double E)
{
    std::transform(m_V.begin(), m_V.end(), m_f.begin(),
            [E, this](double v) {return m_h*m_h*(E - v) / 6 + 1;}
    );

    m_u.front() = 1.0;
    m_u[1] = m_u.front() * (1 + m_h * m_h * ((E - m_V[2])/12 + 7*(m_V[0] - E)/12 - m_h*m_h*(m_V[2] - E)*(m_V[0] - E)/9))
                         / (1 + m_h*m_h*((E - m_V[1])/2 + 2*m_h*m_h*(m_V[1] - E)*(m_V[2] - E)/9));

    for (size_t n = 2; n < m_u.size(); ++n)
    {
        m_u[n] = (m_u[n-1]*(12 - 10*m_f[n-1]) - m_u[n-2]*m_f[n-2])/m_f[n];
    }
}

void Helium::updateDensity()
{
    double E1 = -.55;
    double u1 = 1.0;
    double E2 = -.45;
    double u2 = 1.0;

    while(u1*u2 > 0)
    {
        double delta = (E2 - E1) / 2;

        E1 -= delta;
        shoot(E1);
        u1 = m_u.back();
        if (std::abs(u1/m_u.front()) < epsilon)
        {
            m_E = E1;
            m_u /= norm();
            return;
        }

        E2 += delta;
        shoot(E2);
        u2 = m_u.back();
        if (std::abs(u2/m_u.front()) < epsilon)
        {
            m_E = E2;
            m_u /= norm();
            return;
        }
    }

    while(true)
    {
        double En = (E1 + E2)/2;
        shoot(En);
        if (std::abs(m_u.back()/m_u.front()) < epsilon)
        {
            m_E = En;
            m_u /= norm();
            return;
        }
        else if (u1*m_u.back() > 0)
        {
            u1 = m_u.back();
            E1 = En;
        }
        else{
            E2 = En;
            u2 = m_u.back();
        }
    }
}

std::vector<double> Helium::electrostatic()
{
    if (std::abs(norm() - 1.0) > epsilon*epsilon)
    {
        throw std::runtime_error("State not normalised.");
    }

    for (size_t r = 0; r < m_u.size(); ++r)
    {
        m_f[r] = - m_h*m_u[r]*m_u[r] / (12*(r+1));
    }

    std::vector<double> U(m_u.size(), 0.0);
    U[1] = m_h + (7*m_f[0] + 6*m_f[1] - m_f[2]) / 2;

    for (size_t n = 2; n < U.size(); ++n)
    {
        U[n] = 2*U[n-1] - U[n-2] + m_f[n] + 10*m_f[n-1] + m_f[n-2];
    }

    double k = (1 - U.back()) / U.size();
    for (size_t r = 0; r < U.size(); ++r)
    {
        U[r] += k*(r + 1);
    }

    return U;
}

void Helium::LDA_DFT()
{
    while(true)
    {
        std::vector<double> comparison(m_u);
        std::vector<double> U = electrostatic();
        for (size_t n = 0; n < m_u.size(); ++n)
        {
            m_f[n] = std::pow(m_u[n] / ((n+1)*m_h), 2.0/3.0);
        }
        for (size_t n = 0; n < m_V.size(); ++n)
        {
            m_V[n] = 2 * (U[n] - 1) / ((n + 1) * m_h) - F * m_f[n] - A * std::log(1 + B * m_f[n]);
        }

        updateDensity();
        if (::norm(m_u - comparison, m_h) < epsilon) return;
    }
}
