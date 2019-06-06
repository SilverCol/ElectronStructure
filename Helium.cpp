//
// Created by mitja on 4.6.2019.
//

#include <algorithm>
#include <cmath>
#include <numeric>
#include "Helium.h"


namespace
{
    static const double epsilon = 1e-10;

    double norm(const std::vector<double>& v)
    {
        return std::sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.0));
    }

    void operator/=(std::vector<double>& v, double a)
    {
        for (double& x : v) x /= a;
    }
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
}

void Helium::shoot(double E)
{
    std::transform(m_V.begin(), m_V.end(), m_f.begin(), [E, this](double v){return m_h*m_h*(E-v)/6 + 1;});

    m_u.front() = 1.0;
    m_u[1] = m_u.front()*(12 - 10*m_f.front())/m_f[1];

    for (size_t n = 2; n < m_u.size(); ++n) m_u[n] = (m_u[n-1]*(12 - 10*m_f[n-1]) - m_u[n-2]*m_f[n-2])/m_f[n];
}

void Helium::updateDensity()
{
    double E1 = 0.0;
    shoot(E1);
    double u1 = m_u.back();
    if (std::abs(u1/m_u.front()) < epsilon)
    {
        m_E = E1;
        m_u /= norm(m_u);
        return;
    }

    double E2 = .667868758;
    shoot(E2);
    double u2 = m_u.back();
    if (std::abs(u2/m_u.front()) < epsilon)
    {
        m_E = E2;
        m_u /= norm(m_u);
        return;
    }

    while(u1*u2 > 0)
    {
        E2 *= 2;
        shoot(E2);
        u2 = m_u.back();
        if (std::abs(u2/m_u.front()) < epsilon)
        {
            m_E = E2;
            m_u /= norm(m_u);
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
            m_u /= norm(m_u);
            return;
        }
        else if (u1*m_u.back() > 0)
        {
            u1 = m_u.back();
            E1 = En;
        }
        else E2 = En;
    }
}
