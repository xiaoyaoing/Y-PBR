#pragma once

#include "../Common/math.hpp"
#include <optional>
#include "glm/vec3.hpp"
#include "glm/mat4x4.hpp"
#include "nlohmann/json.hpp"
#include "../Ray/Ray.hpp"
#include "../Ray/Intersection.hpp"
#include "../Common/BoundingBox.hpp"
#include "../Common/util.hpp"
#include "../Lights/Light.hpp"

class AreaLight;

class Bsdf;

class Primitive{
public:
        Primitive(std::shared_ptr<Bsdf> bsdf)
                : bsdf(bsdf) {
//            computeArea();
//            computeBoundingBox();
        };

        virtual ~Primitive() { }

        virtual std::optional<Intersection> intersect(Ray& ray) const = 0;

        virtual vec3 operator()(Float u, Float v) const = 0;

        virtual vec3 normal(const vec3& pos) const = 0;

        virtual void transform(const Transform &T) = 0;

//        virtual const  std::shared_ptr<AreaLight> GetAreaLight() const {
//            return areaLight;
//        }


        // Sample a point on the shape given a reference point |ref| and
        //     return the PDF with respect to solid angle from |ref|.
        virtual Intersection Sample(const Intersection &ref, const vec2 &u,
                                   Float *pdf) const ;

        virtual Intersection Sample(const vec2 &u, Float *pdf) const = 0;

        virtual vec3 interpolatedNormal(const glm::dvec2& uv) const
        {
            return vec3();
        }

        Bounds3 BB() const
        {
            return BB_;
        }

        Float  Area() const {
            return  area;
        }

        Float InvArea() const {
            return inv_area;
        }


        std::shared_ptr<Bsdf> bsdf;

        std::shared_ptr<AreaLight> areaLight;

    int id;
protected:
        virtual void computeArea() = 0;
        virtual void computeBoundingBox() = 0;
        Float area;
        Float inv_area;

        Bounds3 BB_;

    };



