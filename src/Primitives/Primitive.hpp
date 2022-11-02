#pragma once

#include "../Common/math.hpp"
#include "../Common/BoundingBox.hpp"
#include "../Common/Json.hpp"
#include "Common/Frame.hpp"
#include "../Ray/Ray.hpp"
#include "../Ray/Intersection.hpp"
#include "EmbreeUtils.hpp"

#include "../Lights/Light.hpp"


#include <optional>
#include "glm/vec3.hpp"
#include "glm/mat4x4.hpp"
#include "Common/Json.hpp"
#include "Common/Transform.hpp"

class AreaLight;

class BSDF;

class Primitive {
public:
    Primitive(std::shared_ptr < BSDF > bsdf)
            : bsdf(bsdf) {
    };



    virtual ~Primitive( ) {}
    virtual std::optional < Intersection > intersect(Ray & ray) const = 0;
    virtual bool occluded(const Ray & ray) const = 0;
    virtual vec3 operator ()(Float u, Float v) const = 0;
    virtual vec3 normal(const vec3 & pos) const = 0;
    virtual void transform(const mat4 & T) = 0;


    // Sample a point on the shape given a reference point |ref| and
    //     return the PDF with respect to solid angle from |ref|.
    virtual Intersection sample(const Intersection & ref, const vec2 & u,
                                Float * pdf) const;
    virtual Intersection sample(const vec2 & u, Float * pdf) const = 0;
    virtual Float powerToRadianceScale() const {return Constant::INV_PI * inv_area;}
    virtual  Float directPdf(const Intersection & pShape,vec3 ref) const {throw("This function not implemented!");
    }
    virtual Frame setTangentFrame(const Intersection * its) const {
        //todo add bump map
        return Frame(its->Ns);
    }
    virtual vec3 interpolatedNormal(const glm::dvec2 & uv) const {
        return vec3();
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


    RTCGeometry initRTC();

    std::shared_ptr < BSDF > bsdf;
    std::shared_ptr < AreaLight > areaLight;

protected:
    virtual void computeArea( ) = 0;

    virtual void computeBoundingBox( ) = 0;

    Float area;
    Float inv_area;
    Bounds3 BB_;
    RTCGeometry _geom;

    //mat4  toWorld;
};



