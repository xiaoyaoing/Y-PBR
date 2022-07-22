#pragma  once


#include <nlohmann/json.hpp>
#include "../scene.hpp"
#include "../Common/math.hpp"
#include "../Sampler/Sampler.hpp"
#include "../Ray/Ray.hpp"

class Integrator{
public:
   Integrator(nlohmann::json j);

   virtual vec3  integrate(const Ray& ray, const Scene& scene, Sampler& sampler) const =0;
};