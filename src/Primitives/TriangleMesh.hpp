#pragma  once

#include "Primitive.hpp"
#include "TriangleHelper.hpp"
#include "Sampler/Distrib.hpp"
#include "Common/Transform.hpp"

class TriangleMesh : public Primitive {
public:
    TriangleMesh();
    virtual ~ TriangleMesh() = default;
    Intersection sample(const vec2 & u, Float * pdf) const override;
    std::optional<Intersection> intersect(Ray& ray) const override;
    vec3 normal(const vec3& pos) const override;
    void transform(const mat4 & T) override;
    vec3 operator()(Float u, Float v) const override;
    void computeBoundingBox() override;
    void computeArea() override;


    void Load(const Json j,const Scene & scene,const mat4 & transform);
    Bounds3 getTriBounds(int idx);
    Float getTriArea(int idx);
    int BsdfCount() const {return m_bsdfs.size();}
    std::shared_ptr<BSDF> Bsdf(uint32 idx) const {
        return m_bsdfs[idx];
    }

    bool occluded(const Ray & ray) const override;

protected:
    vec2 uvAt(int triID,Float u,Float v) const ;

    std::vector<Vertex> m_vertexs;
    std::vector<TriangleI> m_tris;
    std::vector<std::shared_ptr<BSDF>> m_bsdfs;
    Distribution1D *  m_triSampler ;
    bool useSoomth;

    RTCScene m_scene;
    RTCGeometry m_geometry;
    unsigned _geomId;
};