#include "Primitive.hpp"
#include "Hierancy/BVH.hpp"

#include "IO/CurveIO.hpp"

class Curve : public Primitive{
protected:
    void computeArea( ) override;

    void computeBoundingBox( ) override;

public:
    Curve(const Json & json,Scene & scene);
    void transform(const mat4 & T) override;
    std::optional < Intersection > intersect(Ray & ray) const override;
    bool occluded(const Ray & ray) const override;

    Frame setTangentFrame(const Intersection * its) const override;

    enum CurveModeEnum
    {
        MODE_HALF_CYLINDER,
        MODE_BCSDF_CYLINDER,
        MODE_CYLINDER,
        MODE_RIBBON
    };
private:



    std::vector<uint32> _curveEnds;
    std::vector<vec4> _nodeData;
    std::vector<vec3> _nodeColor;
    std::vector<vec3> _nodeNormals;

    std::vector<uint32> _indices;

    RTCScene m_scene;
    RTCGeometry m_geometry;

    CurveModeEnum _mode;
    Float _subSample;
    Float _curveThickness;
    bool _overrideThickness;

    uint32 _curveCount;
    std::unique_ptr<BVHAccel> _bvh;
};

class CurveI : public  Primitive {
    std::vector<vec4> * _nodeData;
    uint32 id;
    Curve::CurveModeEnum _mode;
    std::optional < Intersection > intersect(Ray & ray) const override;
    bool occluded(const Ray & ray) const override;
    void computeBoundingBox( ) override;

public:
    CurveI(std::vector<vec4> * _nodeData,uint32 id) : _nodeData(_nodeData),id(id){
        computeBoundingBox();
    }
};