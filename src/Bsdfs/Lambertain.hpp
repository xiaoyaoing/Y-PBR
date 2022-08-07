#include "Bsdf.hpp"

class LambertainBsdf : public  Bsdf{

public:

    virtual Spectrum f(const vec3 & wo,const vec3 & wi) const override;

    LambertainBsdf(Spectrum & albedo);
private:

    Spectrum albedo;
};

std::shared_ptr<LambertainBsdf>  CreateLambertainBsdf(nlohmann::json);
