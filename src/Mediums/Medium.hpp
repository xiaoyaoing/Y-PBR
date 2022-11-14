#include "Ray/Ray.hpp"
#include "Colors/Spectrum.hpp"
class  Medium {
public:
    virtual Float TR(const Ray & ray) const =  0;
};

class Homogeneous : public  Medium {
    Float sigmaT = 1.5;
public:
    Float TR(const Ray & ray) const override;
};
