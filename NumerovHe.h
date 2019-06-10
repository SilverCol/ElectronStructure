//
// Created by mitja on 10.6.2019.
//

#ifndef VAJA_VI_1_NUMEROVHE_H
#define VAJA_VI_1_NUMEROVHE_H

#include "Helium.h"

class NumerovHe: public Helium
{
public:
    NumerovHe(const std::vector<double>& u, double R);
private:
    void shoot(double e);
    void updateDensity() override;

};


#endif //VAJA_VI_1_NUMEROVHE_H
