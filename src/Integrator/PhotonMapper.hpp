#include "Integrator.hpp"

class PhotonMapper : public  Integrator{
public:
    void render(const Scene & scene) const override;
    void process(const Scene & scene, Sampler & sampler) override;

    PhotonMapper(const std::shared_ptr < Camera > & camera, int iterations, Float initRadius, int maxBounces,
                 int photonsPerIteration, int writeFrequency);

protected:
    std::shared_ptr<Camera> _camera;
    int iterations;
    Float initRadius;
    int maxBounces;
    int photonsPerIteration;
    int writeFrequency;
};