#pragma once

#include "Common/math.hpp"
#include "Common/BoundingBox.hpp"
#include "Common/Frame.hpp"
#include "Common/Json.hpp"
#include "Common/Transform.hpp"

#include "EmbreeUtils.hpp"

#include "Ray/Ray.hpp"
#include "Ray/Intersection.hpp"
#include "Mediums/Medium.hpp"
#include "Lights/Light.hpp"

#include <optional>


class AreaLight;

class BSDF;

class Primitive {
public:
    Primitive(const Json & json) {
        toWorld = getOptional(json, "transform", getIndentifyTransform());
    }

    Primitive( ) {}

    Primitive(const Bounds3 & bounds, uint32 id) : BB_(bounds), primId(id) {}

    virtual ~Primitive( ) {}

    virtual std::optional < Intersection > intersect(Ray & ray) const {};

    virtual bool occluded(const Ray & ray) const {};

    virtual vec3 operator ()(Float u, Float v) const {};

    virtual vec3 normal(const vec3 & pos) const {};

    virtual void transform(const mat4 & T) {};

    void transform( ) { transform(toWorld); }

    virtual Intersection sample(const vec3 & ref, const vec2 & u,
                                Float * pdf) const;

    virtual Intersection sample(const vec2 & u, Float * pdf) const {};

    virtual Float powerToRadianceScale( ) const { return inv_area; }

    virtual Float directPdf(const Intersection & pShape, vec3 ref) const;

    virtual Frame setTangentFrame(const Intersection * its) const {
        //todo add bump map
        if(its->tangent){
           return Frame(*its->tangent,cross(its->Ns,*its->tangent),its->Ns);
        }
        return Frame(its->Ns);
    }

    const Medium * selectMedium(const Medium * currentMedium, bool geomBack) const  {
        if(inMedium || outMedium)
            return geomBack?inMedium.get():outMedium.get();
        return currentMedium;
    }

    const mat4 & getTransform( ) {
        return toWorld;
    }

    Bounds3 BB( ) const {
        return BB_;
    }

    Float Area( ) const {
        return area;
    }

    Float InvArea( ) const {
        return inv_area;
    }

    virtual  void load(const Json & json, const Scene & scene );
    virtual  bool sameBSDF(BSDF * _bsdf) { return bsdf.get() == _bsdf;}

    RTCGeometry initRTC( );

    std::shared_ptr < BSDF > bsdf;
    std::shared_ptr < BSSRDF > bssrdf;
    std::shared_ptr < AreaLight > areaLight;

    uint32 primId;
protected:
    virtual void computeArea( ) {};

    virtual void computeBoundingBox( ) {};

    Float area;
    Float inv_area;
    Bounds3 BB_;
    RTCGeometry _geom;
    mat4 toWorld;

    std::shared_ptr <Medium > outMedium;
    std::shared_ptr <Medium > inMedium;
    //mat4  toWorld;
};



