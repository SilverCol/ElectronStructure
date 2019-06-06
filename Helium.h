//
// Created by mitja on 4.6.2019.
//

#ifndef VAJA_VI_1_HELIUM_H
#define VAJA_VI_1_HELIUM_H

#include <vector>
#include <cmath>
#include <numeric>

inline void operator/=(std::vector<double>& v, double a)
{
    for (double& x : v) x /= a;
}

inline double norm(const std::vector<double>& v, double h = 1){
    return std::sqrt(h * std::inner_product(v.begin(), v.end(), v.begin(), 0.0));
}

class Helium
{
public:
    Helium(const std::vector<double>& u, double R);
    inline double energy() const {return m_E;}
    inline double norm() const {return ::norm(m_u, m_h);}

private:
    void shoot(double E);
    void updateDensity();
    std::vector<double> electrostatic();

    std::vector<double> m_u;
    std::vector<double> m_V;
    std::vector<double> m_f;
    double m_h;
    double m_E;

};


#endif //VAJA_VI_1_HELIUM_H
