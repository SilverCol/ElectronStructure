//
// Created by mitja on 10.6.2019.
//

#include "NumerovHe.h"

namespace
{
    const double tolerance = 1e-12;
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

bool NumerovHe::nodeless()
{
    return std::find_if(m_u.begin() + 1, m_u.end(),
                [](double u){return u < 0;}) == m_u.end();
}

void NumerovHe::updateDensity()
{
    double E1 = -1.0;
    double E2 = 0.0;

    while(true)
    {
        double En = (E1 + E2)/2;
        shoot(En);
        if (std::abs((E2 - E1)) < tolerance)
        {
            m_e = En;
            m_u /= norm();
            return;
        }
        else if (nodeless())
        {
            E1 = En;
        }
        else{
            E2 = En;
        }
    }
}
