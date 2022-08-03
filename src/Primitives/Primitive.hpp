
#include "../Common/math.hpp"

#pragma once

#include "glm/vec3.hpp"
#include "glm/mat4x4.hpp"
#include "nlohmann/json.hpp"
#include "../Ray/Ray.hpp"
#include "../Ray/Intersection.hpp"
#include "../Common/BoundingBox.hpp"
#include "../Common/util.hpp"
#include "../Lights/Light.hpp"

class Bsdf;

    class Primitive
    {
    public:
        Primitive(std::shared_ptr<Bsdf> bsdf)
                : bsdf(bsdf), area_(0) { };

        virtual ~Primitive() { }

        virtual std::optional<Intersection> intersect(Ray& ray) const = 0;

        virtual vec3 operator()(double u, double v) const = 0;

        virtual vec3 normal(const vec3& pos) const = 0;

        virtual void transform(const Transform &T) = 0;

        virtual const AreaLight *GetAreaLight() const {
            return
        }


        // Sample a point on the shape given a reference point |ref| and
        //     return the PDF with respect to solid angle from |ref|.
        virtual Intersection Sample(const Intersection &ref, const vec2 &u,
                                   Float *pdf) const ;

        virtual Intersection Sample(const vec2 &u, Float *pdf) const = 0;

        virtual vec3 interpolatedNormal(const glm::dvec2& uv) const
        {
            return vec3();
        }

        BoundingBox BB() const
        {
            return BB_;
        }

        double area() const
        {
            return area_;
        }

        std::shared_ptr<Bsdf> bsdf;

    protected:
        virtual void computeArea() = 0;
        virtual void computeBoundingBox() = 0;
        double area_;
        BoundingBox BB_;
        std::shared_ptr<AreaLight> areaLight;
    };

    class Sphere : public Primitive
    {
    public:
        Sphere(double radius, std::shared_ptr<Bsdf> bsdf);

        Intersection Sample(const vec2 & u, Float * pdf) const override;

        virtual std::optional<Intersection> intersect(Ray& ray) const ;

        virtual vec3 operator()(double u, double v) const;

        virtual vec3 normal(const vec3& pos) const;

        virtual void transform(const Transform &T);



    protected:
        virtual void computeArea();
        virtual void computeBoundingBox();

    private:
        vec3 origin;
        Float radius;
    };



    class Quadric : public Primitive
    {
    public:
        Quadric(const nlohmann::json &j, std::shared_ptr<Bsdf> material);

        virtual std::optional < Intersection > intersect(Ray & ray) const;
        virtual vec3 operator()(double u, double v) const;
        virtual vec3 normal(const vec3& pos) const;
        virtual void transform(const Transform &T);

    protected:
        virtual void computeArea();
        virtual void computeBoundingBox() { }

    private:
        glm::dmat4x4 Q; // Quadric matrix
        glm::dmat4x3 G; // Gradient matrix
    };

