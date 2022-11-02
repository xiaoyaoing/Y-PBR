////2022/8/1
//#pragma once
//#include "../Common/math.hpp"
//#include "Primitive.hpp"
//#include "TriangleHelper.hpp"
//
//class Scene;
//class Triangle;
//
//
//struct TriangleMesh{
//    friend  class  Triangle ;
//public:
//    TriangleMesh(std::vector<Vertex> & vertexs){
//        m_vertexs = std::move(vertexs);
//    }
//
//    void Load(const Json j,const Scene & scene,const Transform * transform);
//
//    int BsdfCount() const {return m_bsdfs.size();}
//
//    std::shared_ptr<BSDF> Bsdf(uint32 idx) const {
//        return m_bsdfs[idx];
//    }
//private:
//    std::vector<Vertex> m_vertexs;
//    std::vector<std::shared_ptr<BSDF>> m_bsdfs;
//    bool useSoomth;
//};
//
//
//class Triangle : public Primitive
//{
//public:
//    Triangle(const std::shared_ptr<TriangleMesh> &mesh,
//             const std::shared_ptr<BSDF> & bsdf,
//             uint32 * v_indices)
//             : Primitive(bsdf), mesh(mesh)
//             {
//             m_v=v_indices;
//             computeArea();
//             computeBoundingBox();
//    }
//
//    virtual std::optional < Intersection > intersect(Ray & ray) const;
//
//    virtual vec3 operator()(Float u, Float v) const;
//
//    virtual vec3 normal(const vec3& pos) const;
//
//    virtual vec3 interpolatedNormal(const glm::dvec2& uv) const;
//
//    virtual void transform(const Transform &T);
//
//
//
//    Intersection Sample(const vec2 & u, Float * pdf) const override;
//
//    vec3 normal(const vec2 & uv) const ;
//    vec3 selectNs(const vec2 & uv,const vec3 & Ng) const ;
//protected:
//    virtual void computeArea();
//    virtual void computeBoundingBox();
//
//    uint32 * m_v;
////        const std::unique_ptr<glm::dmat3> N; // vertex normals
//    std::shared_ptr <TriangleMesh> mesh;
//
//};
//
//std::shared_ptr<TriangleMesh> CreateTriangleMesh(
////        const Transform *o2w, const Transform *w2o, bool reverseOrientation,
//        const  Transform * transform,
//      //  const  std::vector<std::vector<size_t>> * vertexIndices,
//        const std::vector<vec3>*p,
//        const std::vector<vec3> *s, const std::vector<vec3> *n, const std::vector<vec2> *uv,
//     //   std::shared_ptr<Bsdf> bsdf,
//        const std::vector<size_t> * faceIndices = nullptr);
//
//
//std::vector<std::shared_ptr<Primitive>> getTrianglesFromMesh(
//                                        const std::shared_ptr<TriangleMesh> mesh,
//                                        const std::vector<TriangleI> & tris);
//
