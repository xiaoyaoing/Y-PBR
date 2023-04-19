#include "Integrator.hpp"

class PhotonMapper : public  Integrator{
public:
    void render(const Scene & scene)  override;
    void process(const Scene & scene, Sampler & sampler) override;

    PhotonMapper(const std::shared_ptr < Camera > & camera, const Json &);

protected:
    std::shared_ptr<Camera> _camera;
    int iterations;
    Float initRadius;
    int photonsPerIteration;
    int writeFrequency;
};