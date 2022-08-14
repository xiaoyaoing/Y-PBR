//2022/7/17

#include "Triangle.hpp"
#include "spdlog/spdlog.h"

TriangleMesh::TriangleMesh(const Transform * T,
                           //const std::vector < std::vector < size_t>> * VertexIndices,
                           const std::vector < vec3 > * P, const std::vector < vec3 > * S,
                           const std::vector < vec3 > * N, const std::vector < vec2 > * UV)
                           {
    int nVertices = P->size();
    p.reset(new vec3[nVertices]);

    for (int i = 0; i < nVertices; ++i) p[i] =T?(*T)(P->at(i)):P->at(i);
    //

    if(S){

    }
    //
    if(N){

    }
    //
    if(UV){

    }


}

std::optional < Intersection > Triangle::intersect(Ray & ray) const {
    /* Begin calculating determinant - also used to calculate U parameter */
    Intersection intersection;

    const vec3  &p0 = mesh->p[v[0]];
    const vec3  &p1 = mesh->p[v[1]];
    const vec3  &p2 = mesh->p[v[2]];

    const vec3 e1  = p1-p0;
    const vec3 e2  = p2-p0;

    Float u,v,t;
    vec3 pvec = cross(ray.d,e2);

    /* If determinant is near zero, ray lies in plane of triangle */
    float det = dot(e1,pvec);

    if (det > -1e-8f && det < 1e-8f)
        return  std::nullopt;
    float inv_det = 1.0f / det;

    /* Calculate distance from v[0] to ray origin */
    vec3 tvec = ray.o - p0;


    /* Calculate U parameter and test bounds */
    u =dot(pvec,tvec) * inv_det;
    if (u < 0.0 || u > 1.0)
        return std::nullopt;

    /* Prepare to test V parameter */
    vec3 qvec = cross(tvec,e1);

    /* Calculate V parameter and test bounds */
    v = dot(ray.d,qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0)
        return std::nullopt;

    /* Ray intersects triangle -> compute t */
    t = dot(e2,qvec) * inv_det;

    if(t >= ray.nearT && t <= ray.farT){
        ray.farT=t;

        intersection.setNormal(-normalize(cross(e1, e2)));

        intersection.p=ray(t);
        intersection.primitive=this;
        intersection.bsdf=bsdf.get();
        return {intersection};
    }
    return std::nullopt;
}

vec3 Triangle::operator ()(double u, double v) const {
    return vec3(0);
}

vec3 Triangle::normal(const vec3 & pos) const {
    const vec3  &p0 = mesh->p[v[0]];
    const vec3  &p1 = mesh->p[v[1]];
    const vec3  &p2 = mesh->p[v[2]];

    const vec3 e1  = p1-p0;
    const vec3 e2  = p2-p0;

    return normalize(cross(e1,e2));
}

void Triangle::transform(const Transform & T) {
     spdlog::error("Should not transform single Triangle");
//    v0 = T.matrix * vec4(v0, 1.0);
//    v1 = T.matrix * vec4(v1, 1.0);
//    v2 = T.matrix * vec4(v2, 1.0);
//
//    e1 = v1 - v0;
//    e2 = v2 - v0;
//    normal_ = normalize(cross(e1, e2));
}


vec3 Triangle::interpolatedNormal(const glm::dvec2 & uv) const {
    return Primitive::interpolatedNormal(uv);
}

void Triangle::computeArea( ) {
    const vec3 &p0 = mesh->p[v[0]];
    const vec3 &p1 = mesh->p[v[1]];
    const vec3 &p2 = mesh->p[v[2]];
    area=0.5 * glm::length(cross(p1 - p0, p2 - p0));
    inv_area=1/area;
}

void Triangle::computeBoundingBox( ) {

}

Intersection Triangle::Sample(const vec2 & u, Float * pdf) const {
    Intersection it;

    const vec3  &p0 = mesh->p[v[0]];
    const vec3  &p1 = mesh->p[v[1]];
    const vec3  &p2 = mesh->p[v[2]];

    vec3 point = p0*u.x + p1*u.y +(1-u.x-u.y)*p2;

    it.p=point;
    it.setNormal(normalize((cross(p1 - p0, p2 - p0))));
    if(mesh->n){
        //todo support mesh n
    }
    *pdf=inv_area;

    return it;
}



std::shared_ptr<TriangleMesh> CreateTriangleMesh(
//        const Transform *o2w, const Transform *w2o, bool reverseOrientation,
        const Transform * transform,
       // const  std::vector<std::vector<size_t>> * vertexIndices,
        const std::vector<vec3>*p,
        const std::vector<vec3> *s, const std::vector<vec3> *n, const std::vector<vec2> *uv,
        //std::shared_ptr<Bsdf> bsdf,
        const std::vector<size_t> * faceIndices){

    auto mesh =std::make_shared<TriangleMesh>(transform,p,s,n,uv);

    std::vector<std::shared_ptr<Primitive>> primitives;

    return mesh;



}

std::vector<std::shared_ptr<Primitive>> getTrianglesFromMesh(
                                        const std::shared_ptr<TriangleMesh> mesh,
                                        const std::vector<std::vector<size_t>> v_indces,
                                        const std::shared_ptr<Bsdf> bsdf
                                        ){

    std::vector<std::shared_ptr<Primitive>> triangles;

    for(int i=0;i<v_indces.size();i++){
        std::shared_ptr<Primitive> triangle = std::make_shared<Triangle>(mesh,bsdf,v_indces[i]);
        triangles.push_back(triangle);
    }

    return triangles;
}





