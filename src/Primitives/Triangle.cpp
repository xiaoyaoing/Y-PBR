////2022/7/17
//
//#include "Triangle.hpp"
//#include <spdlog/spdlog.h>
//#include "scene.hpp"
//
//
//void TriangleMesh::Load(const Json j, const Scene & scene,const Transform * transform) {
//    auto & bsdf_json = j["bsdf"];
//
//    if(bsdf_json.is_array()){
//        for(const std::string & bsdf_str : bsdf_json){
//            m_bsdfs.push_back(scene.fetchBSDF(bsdf_str));
//        }
//    }
//    else {
//        m_bsdfs.push_back(scene.fetchBSDF(bsdf_json));
//    }
//
//    useSoomth = getOptional(j,"use_smooth",true);
//
//    if(transform == nullptr)
//        return ;
//
//    mat4  transformNormalMat = transform->TransformNormalMat();
//
//    for(Vertex & vertex:m_vertexs){
//        vertex.pos() = transform ->operator *(vertex.pos());
//        vertex.normal() = mult(transformNormalMat,vec4(vertex.normal(),0));
//    }
//}
//
//std::optional < Intersection > Triangle::intersect(Ray & ray) const {
//    /* Begin calculating determinant - also used to calculate U parameter */
//    Intersection intersection;
//
//
//
//    const vec3 & p0 = mesh->m_vertexs[m_v[0]].pos();
//    const vec3 & p1 = mesh->m_vertexs[m_v[1]].pos();
//    const vec3 & p2 = mesh->m_vertexs[m_v[2]].pos();
//
//    const vec3 e1 = p1 - p0;
//    const vec3 e2 = p2 - p0;
//
//    Float u, v, t;
//    vec3 pvec = cross(ray.d, e2);
//
//    /* If determinant is near zero, ray lies in plane of triangle */
//    float det = dot(e1, pvec);
//
//    if ( det > - 1e-8f && det < 1e-8f )
//        return std::nullopt;
//    float inv_det = 1.0f / det;
//
//    /* Calculate distance from v[0] to ray origin */
//    vec3 tvec = ray.o - p0;
//
//
//    /* Calculate U parameter and test bounds */
//    u = dot(pvec, tvec) * inv_det;
//    if ( u < 0.0 || u > 1.0 )
//        return std::nullopt;
//
//    /* Prepare to test V parameter */
//    vec3 qvec = cross(tvec, e1);
//
//    /* Calculate V parameter and test bounds */
//    v = dot(ray.d, qvec) * inv_det;
//    if ( v < 0.0 || u + v > 1.0 )
//        return std::nullopt;
//
//    /* Ray intersects triangle -> compute t */
//    t = dot(e2, qvec) * inv_det;
//
//    if ( t >= ray.nearT && t <= ray.farT ) {
//        ray.farT = t;
//        intersection.Ng=- normalize(cross(e1, e2));
//        intersection.Ns = selectNs(vec2(u,v),intersection.Ng);
//        intersection.p = ray(t);
//        intersection.primitive = this;
//        intersection.bsdf = bsdf.get();
//        return {intersection};
//    }
//    return std::nullopt;
//}
//
//vec3 Triangle::operator ()(Float u, Float v) const {
//    return vec3(0);
//}
//
//vec3 Triangle::normal(const vec3 & pos) const {
//    const vec3 & p0 = mesh->m_vertexs[m_v[0]].pos();
//    const vec3 & p1 = mesh->m_vertexs[m_v[1]].pos();
//    const vec3 & p2 = mesh->m_vertexs[m_v[2]].pos();
//
//    const vec3 e1 = p1 - p0;
//    const vec3 e2 = p2 - p0;
//
//    return normalize(cross(e1, e2));
//}
//
//void Triangle::transform(const Transform & T) {
//    spdlog::error("Should not transform single Triangle");
////    v0 = T.matrix * vec4(v0, 1.0);
////    v1 = T.matrix * vec4(v1, 1.0);
////    v2 = T.matrix * vec4(v2, 1.0);
////
////    e1 = v1 - v0;
////    e2 = v2 - v0;
////    normal_ = normalize(cross(e1, e2));
//}
//
//
//vec3 Triangle::interpolatedNormal(const glm::dvec2 & uv) const {
//    return Primitive::interpolatedNormal(uv);
//}
//
//void Triangle::computeArea( ) {
//    const vec3 & p0 = mesh->m_vertexs[m_v[0]].pos();
//    const vec3 & p1 = mesh->m_vertexs[m_v[1]].pos();
//    const vec3 & p2 = mesh->m_vertexs[m_v[2]].pos();
//    area = 0.5 * glm::length(cross(p1 - p0, p2 - p0));
//    inv_area = 1 / area;
//}
//
//void Triangle::computeBoundingBox( ) {
//    const vec3 & p0 = mesh->m_vertexs[m_v[0]].pos();
//    const vec3 & p1 = mesh->m_vertexs[m_v[1]].pos();
//    const vec3 & p2 = mesh->m_vertexs[m_v[2]].pos();
//
//    vec3 pMin = min(p0, min(p1, p2));
//    vec3 pMax = max(p0, max(p1, p2));
//
//    BB_ = Bounds3(pMin, pMax);
//}
//
//
//Intersection Triangle::Sample(const vec2 & u, Float * pdf) const {
//    Intersection it;
//
//    const vec3 & p0 = mesh->m_vertexs[m_v[0]].pos();
//    const vec3 & p1 = mesh->m_vertexs[m_v[1]].pos();
//    const vec3 & p2 = mesh->m_vertexs[m_v[2]].pos();
//
//    vec3 point = p0 * u.x + p1 * u.y + ( 1 - u.x - u.y ) * p2;
//
//    it.p = point;
//    it.Ng=normalize(( cross(p1 - p0, p2 - p0) ));
//
//    it.Ns = selectNs(u,it.Ng);
//
//    * pdf = inv_area;
//
//    return it;
//}
//
//
//vec3 Triangle::normal(const vec2 & uv) const  {
//    return  uv.x * mesh->m_vertexs[m_v[0]].normal()  + uv.y * mesh->m_vertexs[m_v[1]].normal()
//            + (1-uv.x-uv.y) * mesh->m_vertexs[m_v[2]].normal();
//}
//
//vec3 Triangle::selectNs(const vec2 & uv, const vec3 & Ng) const {
//    if(mesh->useSoomth){
//        return normal(uv);
//    }
//    return Ng;
//}
//
//
//
////std::shared_ptr < TriangleMesh > CreateTriangleMesh(
//////        const Transform *o2w, const Transform *w2o, bool reverseOrientation,
////        const Transform * transform,
////        // const  std::vector<std::vector<size_t>> * vertexIndices,
////        const std::vector < vec3 > * p,
////        const std::vector < vec3 > * s, const std::vector < vec3 > * n, const std::vector < vec2 > * uv,
////        //std::shared_ptr<Bsdf> bsdf,
////        const std::vector < size_t > * faceIndices) {
////
////        auto mesh = std::make_shared < TriangleMesh >(transform, p, s, n, uv);
////
////        return mesh;
////
////
////}
//
//std::vector < std::shared_ptr < Primitive>> getTrianglesFromMesh(
//        const std::shared_ptr < TriangleMesh > mesh,
//        const std::vector<TriangleI> &  tris
//) {
//
//    std::vector < std::shared_ptr < Primitive>> triangles;
//
//    for ( int i = 0 ; i < tris.size() ; i ++ ) {
//        uint32  indices[3]  = {tris[i].v0,tris[i].v1,tris[i].v2};
//        std::shared_ptr < Primitive > triangle = std::make_shared < Triangle >
//                            (mesh,mesh->Bsdf(clamp(tris[i].material,0,mesh->BsdfCount()-1)),indices);
//        triangles.push_back(triangle);
//    }
//
//    return triangles;
//}
//
//
//
//
//
