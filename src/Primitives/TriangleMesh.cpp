#include "TriangleMesh.hpp"
#include "IO/MeshIO.hpp"
#include "Sampler/Distrib.hpp"
#include "scene.hpp"
#include "Common/Transform.hpp"

TriangleMesh::TriangleMesh( ) : Primitive(nullptr) {}


void TriangleMesh::Load(const Json & json, const Scene & scene) {
    MeshIO::LoadMeshFromFile(json["file"].get < std::string >(), m_vertexs, m_tris);
    mat4 transformMatrix = getOptional(json, "transform", getIndentifyTransform());

    Json bsdf_json = json["bsdf"];
    if ( bsdf_json.is_array() ) {
        for ( const std::string & bsdf_str: bsdf_json ) {
            m_bsdfs.push_back(scene.fetchBSDF(bsdf_str));
        }
    } else {
        m_bsdfs.push_back(scene.fetchBSDF(bsdf_json));
    }

    useSoomth = getOptional(json, "use_smooth", true);
    // useSoomth = false;
    if ( transformMatrix!=mat4() ) {
        mat4 transformNormalMat = getTransformNormalMat(transformMatrix);
        std::string s = Mat4ToStr(transformNormalMat);
        for ( Vertex & vertex: m_vertexs ) {
            vertex.pos() = transformPoint(transformMatrix, vertex.pos());
            vertex.normal() = transformVector(transformNormalMat, vertex.normal());
        }
    }


    computeBoundingBox();
    computeArea();

    m_scene = rtcNewScene(EmbreeUtils::getDevice());

    m_geometry = rtcNewGeometry(EmbreeUtils::getDevice(), RTC_GEOMETRY_TYPE_TRIANGLE);
    // set vertices
    float * vb = (float *) rtcSetNewGeometryBuffer(m_geometry, RTC_BUFFER_TYPE_VERTEX, 0,
                                                   RTC_FORMAT_FLOAT3,
                                                   3 * sizeof(float), m_vertexs.size());
    for ( size_t i = 0 ; i < m_vertexs.size() ; i ++ ) {
        vb[3 * i] = m_vertexs[i].pos().x;
        vb[3 * i + 1] = m_vertexs[i].pos().y;
        vb[3 * i + 2] = m_vertexs[i].pos().z;
    }

    // set indices
    unsigned * ib = (unsigned *) rtcSetNewGeometryBuffer(m_geometry, RTC_BUFFER_TYPE_INDEX, 0,
                                                         RTC_FORMAT_UINT3,
                                                         3 * sizeof(unsigned), m_tris.size());
    for ( size_t i = 0 ; i < m_tris.size() ; i ++ ) {
        ib[i * 3] = m_tris[i].v0;
        ib[i * 3 + 1] = m_tris[i].v1;
        ib[i * 3 + 2] = m_tris[i].v2;
    }


    for ( TriangleI & triangleI: m_tris ) triangleI.material = clamp(triangleI.material, 0, m_bsdfs.size() - 1);

    rtcCommitGeometry(m_geometry);
    rtcAttachGeometry(m_scene, m_geometry);
    // rtcReleaseGeometry(geom);
    rtcCommitScene(m_scene);
}


void TriangleMesh::computeBoundingBox( ) {
    Bounds3 bb;
    for ( int i = 0 ; i < m_tris.size() ; i ++ ) {
        bb = Union(bb, getTriBounds(i));
    }
    BB_ = bb;
}

void TriangleMesh::computeArea( ) {
    area = 0;
    Float * areaDis = new Float[m_tris.size()];
    for ( int i = 0 ; i < m_tris.size() ; i ++ ) {
        Float curArea = getTriArea(i);
        area += curArea;
        areaDis[i] = curArea;
    }
    m_triSampler = new Distribution1D(areaDis, m_tris.size());

}

Bounds3 TriangleMesh::getTriBounds(int idx) {
    auto & m_v = m_tris[idx].vs;
    const vec3 & p0 = m_vertexs[m_v[0]].pos();
    const vec3 & p1 = m_vertexs[m_v[1]].pos();
    const vec3 & p2 = m_vertexs[m_v[2]].pos();

    vec3 pMin = min(p0, min(p1, p2));
    vec3 pMax = max(p0, max(p1, p2));

    return Bounds3(pMin, pMax);
}

Float TriangleMesh::getTriArea(int idx) {
    auto & m_v = m_tris[idx].vs;
    const vec3 & p0 = m_vertexs[m_v[0]].pos();
    const vec3 & p1 = m_vertexs[m_v[1]].pos();
    const vec3 & p2 = m_vertexs[m_v[2]].pos();
    return 0.5 * glm::length(cross(p1 - p0, p2 - p0));
}

std::optional < Intersection > TriangleMesh::intersect(Ray & ray) const {
    RTCRayHit rayHit;
    EmbreeUtils::convertRay(& ray, & rayHit);

    RTCIntersectContext context;
    rtcInitIntersectContext(& context);

    rtcIntersect1(m_scene, & context, & rayHit);
    if ( rayHit.hit.geomID != RTC_INVALID_GEOMETRY_ID ) {
        ray.farT = rayHit.ray.tfar;

        Intersection its;
        its.p = ray.operator ()(ray.farT);
        const TriangleI tri = m_tris[rayHit.hit.primID];
        its.bsdf = Bsdf(tri.material).get();
        its.uv = uvAt(rayHit.hit.primID, rayHit.hit.u, rayHit.hit.v);


        const vec3 & p0 = m_vertexs[tri.v0].pos();
        const vec3 & p1 = m_vertexs[tri.v1].pos();
        const vec3 & p2 = m_vertexs[tri.v2].pos();

        its.Ns = its.Ng = normalize(cross(p1 - p0, p2 - p0));
        if ( useSoomth ) {
            const vec3 & n1 = m_vertexs[tri.v0].normal();
            const vec3 & n2 = m_vertexs[tri.v1].normal();
            const vec3 & n3 = m_vertexs[tri.v2].normal();
            its.Ns = normalize(interpolate3(n2, n3, n1, vec2(rayHit.hit.u, rayHit.hit.v)));
            //  its.Ns = vec3(rayHit.hit.u, rayHit.hit.v,0);
            // its.Ns = normalize(n1);
            // its.Ns = vec3(Float(tri.v0)/100000,(float)tri.v1 / 100000,(float)tri.v2 / 100000);
        }
        its.primitive = this;
        return {its};
    } else {
        return std::nullopt;
    }
}

bool TriangleMesh::occluded(const Ray & ray) const {
    RTCRay rtcRay;
    RTCIntersectContext context;
    rtcInitIntersectContext(& context);
    EmbreeUtils::convertRay(& ray, & rtcRay);
    rtcOccluded1(m_scene, & context, & rtcRay);
    if ( rtcRay.tfar != - std::numeric_limits < Float >::infinity() )
        return false;
    return true;
}

vec3 TriangleMesh::normal(const vec3 & pos) const {
    throw ( "not implmented" );
}

vec3 TriangleMesh::operator ()(Float u, Float v) const {
    throw ( "not implmented" );
}

void TriangleMesh::transform(const mat4 & T) {
    throw ( "not implmented" );
}

Intersection TriangleMesh::sample(const vec2 & u, Float * pdf) const {
    return Intersection();
}

vec2 TriangleMesh::uvAt(int triID, Float u, Float v) const {
    const TriangleI & t = m_tris[triID];
    vec2 uv0 = m_vertexs[t.v0].uv();
    vec2 uv1 = m_vertexs[t.v1].uv();
    vec2 uv2 = m_vertexs[t.v2].uv();
    return ( 1.0f - u - v ) * uv0 + u * uv1 + v * uv2;
}



