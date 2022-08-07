#include "Lambertain.hpp"
#include "../Common/util.hpp"

Spectrum LambertainBsdf::f(const vec3 & wo,const vec3 & wi) const {
    return albedo * Constant::INV_PI;
}

LambertainBsdf::LambertainBsdf(Spectrum & albedo) :albedo(albedo){

}


std::shared_ptr<LambertainBsdf>  CreateLambertainBsdf(nlohmann::json j){
    Spectrum albedo =getOptional(j, "scale", Spectrum(1.0));
    return std::make_shared<LambertainBsdf>(albedo);
}

