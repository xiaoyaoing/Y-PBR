//2022/8/1
#pragma once
#include "../Common/math.hpp"
#include "Primitive.hpp"

struct TriangleMesh{
   // std::vector<std::vector<size_t>>  vertexIndices;
    std::unique_ptr<vec3[]> p;
    std::unique_ptr<vec3[]> n;
    std::unique_ptr<vec3[]> s;
    std::unique_ptr<vec3[]> uv;
   // std::shared_ptr<Bsdf>   bsdf;

    TriangleMesh(  const Transform * transform,
     //              const  std::vector<std::vector<size_t>> * VertexIndices,
                   const std::vector<vec3>* P, const std::vector<vec3> * S,
                   const std::vector<vec3>* N, const std::vector<vec2> * UV
                //   std::shared_ptr<Bsdf> bsdf
                   );



};

class Triangle : public Primitive
{
public:
    Triangle(const std::shared_ptr<TriangleMesh> &mesh,
             const std::shared_ptr<Bsdf> & bsdf,
             const std::vector<size_t> & v_indices)
             : Primitive(bsdf), mesh(mesh)
             {
             v=v_indices;
             computeArea();
             computeBoundingBox();
    }

    virtual std::optional < Intersection > intersect(Ray & ray) const;

    virtual vec3 operator()(double u, double v) const;

    virtual vec3 normal(const vec3& pos) const;

    virtual vec3 interpolatedNormal(const glm::dvec2& uv) const;

    virtual void transform(const Transform &T);

    vec3 normal() const;

    Intersection Sample(const vec2 & u, Float * pdf) const override;


protected:
    virtual void computeArea();
    virtual void computeBoundingBox();

     std::vector<size_t> v;
//        const std::unique_ptr<glm::dmat3> N; // vertex normals
    std::shared_ptr <TriangleMesh> mesh;

};

std::shared_ptr<TriangleMesh> CreateTriangleMesh(
//        const Transform *o2w, const Transform *w2o, bool reverseOrientation,
        const  Transform * transform,
      //  const  std::vector<std::vector<size_t>> * vertexIndices,
        const std::vector<vec3>*p,
        const std::vector<vec3> *s, const std::vector<vec3> *n, const std::vector<vec2> *uv,
     //   std::shared_ptr<Bsdf> bsdf,
        const std::vector<size_t> * faceIndices = nullptr);


std::vector<std::shared_ptr<Primitive>> getTrianglesFromMesh(
                                        const std::shared_ptr<TriangleMesh> mesh,
                                        const std::vector<std::vector<size_t>> v_indces,
                                        const std::shared_ptr<Bsdf> bsdf);

