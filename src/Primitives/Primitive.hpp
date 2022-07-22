
#include "../Common/math.hpp"

#pragma once

#include "glm/vec3.hpp"
#include "glm/mat4x4.hpp"
#include "nlohmann/json.hpp"

#include "../Ray/Ray.hpp"
#include "../Ray/Intersection.hpp"
#include "../Common/BoundingBox.hpp"
#include "../Common/util.hpp"

class Bsdf;

    class Primitive
    {
    public:
        Primitive(std::shared_ptr<Bsdf> bsdf)
                : bsdf(bsdf), area_(0) { };

        virtual ~Primitive() { }

        virtual bool intersect(Ray& ray, Intersection& intersection) const = 0;
        virtual vec3 operator()(double u, double v) const = 0;
        virtual vec3 normal(const vec3& pos) const = 0;
        virtual void transform(const Transform &T) = 0;

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
    };

    class Sphere : public Primitive
    {
    public:
        Sphere(double radius, std::shared_ptr<Bsdf> bsdf);

        virtual bool intersect(Ray& ray, Intersection& intersection) const;
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

    class Triangle : public Primitive
    {
    public:
        Triangle(const vec3& v0, const vec3& v1, const vec3& v2, std::shared_ptr<Bsdf> Bsdf);

        Triangle(const vec3& v0, const vec3& v1, const vec3& v2,
                 const vec3& n0, const vec3& n1, const vec3& n2, std::shared_ptr<Bsdf> Bsdf);

        virtual bool intersect( Ray& ray, Intersection& intersection) const;
        virtual vec3 operator()(double u, double v) const;
        virtual vec3 normal(const vec3& pos) const;
        virtual vec3 interpolatedNormal(const glm::dvec2& uv) const;
        virtual void transform(const Transform &T);

        vec3 normal() const;

    protected:
        virtual void computeArea();
        virtual void computeBoundingBox();

        vec3 v0, v1, v2;
        const std::unique_ptr<glm::dmat3> N; // vertex normals

        // Pre-computed edges and normal
        vec3 E1, E2, normal_;
    };

    class Quadric : public Primitive
    {
    public:
        Quadric(const nlohmann::json &j, std::shared_ptr<Bsdf> material);

        virtual bool intersect( Ray& ray, Intersection& intersection) const;
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

