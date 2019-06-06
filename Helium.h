//
// Created by mitja on 4.6.2019.
//

#ifndef VAJA_VI_1_HELIUM_H
#define VAJA_VI_1_HELIUM_H

#include <vector>

class Helium
{
public:
    Helium(const std::vector<double>& u, double R);
    inline double energy() const {return m_E;}

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
