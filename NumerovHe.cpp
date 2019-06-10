//
// Created by mitja on 10.6.2019.
//

#include "NumerovHe.h"

namespace
{
    const double tolerance = 1e-4;
}

NumerovHe::NumerovHe(const std::vector<double>& u, double R) :
Helium(u, R)
{
    LDA_DFT();
}

void NumerovHe::shoot(double e)
{
    std::transform(m_V.begin(), m_V.end(), m_f.begin(),
                   [this, e](double v) {return m_h*m_h*(e - v) / 6 + 1;}
    );

    m_u.front() = 0.0;
    m_u[1] = 1.0;

    for (size_t n = 2; n < m_u.size(); ++n)
    {
        m_u[n] = (m_u[n-1]*(12 - 10*m_f[n-1]) - m_u[n-2]*m_f[n-2])/m_f[n];
    }
}

void NumerovHe::updateDensity()
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
        if (std::abs(u1/m_u.front()) < tolerance)
        {
            m_e = E1;
            m_u /= norm();
            return;
        }

        E2 += delta;
        shoot(E2);
        u2 = m_u.back();
        if (std::abs(u2/m_u.front()) < tolerance)
        {
            m_e = E2;
            m_u /= norm();
            return;
        }
    }

    while(true)
    {
        double En = (E1 + E2)/2;
        shoot(En);
        if (std::abs((E2 - E1) / En) < tolerance)
        {
            m_e = En;
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
