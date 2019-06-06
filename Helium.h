//
// Created by mitja on 4.6.2019.
//

#ifndef VAJA_VI_1_HELIUM_H
#define VAJA_VI_1_HELIUM_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

inline void operator/=(std::vector<double>& v, double a)
{
    for (double& x : v) x /= a;
}

inline std::vector<double> operator-(const std::vector<double>& v1, const std::vector<double>& v2)
{
    std::vector<double> w(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), w.begin(), std::minus<>());
    return w;
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
    std::vector<double> electrostatic();

private:
    void shoot(double E);
    void updateDensity();
    void LDA_DFT();

    std::vector<double> m_u;
    std::vector<double> m_V;
    std::vector<double> m_f;
    double m_h;
    double m_E;

};


#endif //VAJA_VI_1_HELIUM_H
