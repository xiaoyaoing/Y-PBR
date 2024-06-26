#include "TriangleMesh.hpp"
#include "IO/MeshIO.hpp"
#include "Sampler/Distrib.hpp"
#include "scene.hpp"
#include "Common/Transform.hpp"

#include <unordered_map>

void TriangleMesh::load(const Json& json, const Scene& scene) {
    Primitive::load(json, scene);
    loadResources(json, scene);
    useSoomth             = getOptional(json, "smooth", false);
    bool recomputeNormals = getOptional(json, "recompute_normals", false);
    if (useSoomth && recomputeNormals)
        recomputeSmoothNormals();

    computeBoundingBox();
    computeArea();
    buildRTC();
}

void TriangleMesh::computeBoundingBox() {
    Bounds3 bb;
    for (int i = 0; i < m_tris.size(); i++) {
        bb = Union(bb, getTriBounds(i));
    }

    BB_ = bb;
}

void TriangleMesh::computeArea() {
    area           = 0;
    Float* areaDis = new Float[m_tris.size()];
    for (int i = 0; i < m_tris.size(); i++) {
        Float curArea = getTriArea(i);
        area += curArea;
        areaDis[i] = curArea;
    }
    m_triSampler = new Distribution1D(areaDis, m_tris.size());
}

Bounds3 TriangleMesh::getTriBounds(int idx) {
    auto&       m_v = m_tris[idx].vs;
    const vec3& p0  = m_vertexs[m_v[0]].pos();
    const vec3& p1  = m_vertexs[m_v[1]].pos();
    const vec3& p2  = m_vertexs[m_v[2]].pos();

    vec3 pMin = min(p0, min(p1, p2));
    vec3 pMax = max(p0, max(p1, p2));

    return Bounds3(pMin, pMax);
}

Float TriangleMesh::getTriArea(int idx) {
    auto&       m_v = m_tris[idx].vs;
    const vec3& p0  = m_vertexs[m_v[0]].pos();
    const vec3& p1  = m_vertexs[m_v[1]].pos();
    const vec3& p2  = m_vertexs[m_v[2]].pos();
    return 0.5 * glm::length(cross(p1 - p0, p2 - p0));
}

std::optional<Intersection> TriangleMesh::intersect(Ray& ray) const {
    RTCRayHit rayHit;
    EmbreeUtils::convertRay(&ray, &rayHit);

    RTCIntersectContext context;
    rtcInitIntersectContext(&context);

    rtcIntersect1(m_scene, &context, &rayHit);
    if (rayHit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
        ray.farT = rayHit.ray.tfar;

        Intersection    its;
        its.p               = ray.    operator()(ray.farT);
        const TriangleI tri = m_tris[rayHit.hit.primID];
        its.bsdf            = Bsdf(tri.material).get();
        its.bssrdf          = bssrdf.get();
        its.uv              = uvAt(rayHit.hit.primID, rayHit.hit.u, rayHit.hit.v);

        const vec3& p0 = m_vertexs[tri.v0].pos();
        const vec3& p1 = m_vertexs[tri.v1].pos();
        const vec3& p2 = m_vertexs[tri.v2].pos();

        its.Ns = its.Ng = normalize(cross(p1 - p0, p2 - p0));
        // 定义随机数分布，范围为 [0, 99]
        if (useSoomth) {
            const vec3& n1 = m_vertexs[tri.v0].normal();
            const vec3& n2 = m_vertexs[tri.v1].normal();
            const vec3& n3 = m_vertexs[tri.v2].normal();
            its.Ns         = normalize(interpolate3(n2, n3, n1, vec2(rayHit.hit.u, rayHit.hit.v)));
        }
        if (hasNan(its.Ns)) {
        }
        its.primitive = this;
        return {its};
    } else {
        return std::nullopt;
    }
}

bool TriangleMesh::occluded(const Ray& ray) const {
    RTCRay              rtcRay;
    RTCIntersectContext context;
    rtcInitIntersectContext(&context);
    EmbreeUtils::convertRay(&ray, &rtcRay);
    rtcOccluded1(m_scene, &context, &rtcRay);
    if (rtcRay.tfar != -std::numeric_limits<Float>::infinity())
        return false;
    return true;
}

vec3 TriangleMesh::normal(const vec3& pos) const {
    throw("not implmented");
}

void TriangleMesh::transform(const mat4& T) {
    //throw ( "not implmented" );
}

Intersection TriangleMesh::sample(const vec2& u, Float* pdf) const {
    return Intersection();
}

vec2 TriangleMesh::uvAt(int triID, Float u, Float v) const {
    const TriangleI& t   = m_tris[triID];
    vec2             uv0 = m_vertexs[t.v0].uv();
    vec2             uv1 = m_vertexs[t.v1].uv();
    vec2             uv2 = m_vertexs[t.v2].uv();
    return (1.0f - u - v) * uv0 + u * uv1 + v * uv2;
}

void TriangleMesh::recomputeSmoothNormals() {
    static const float SplitLimit = std::cos(Constant::PI * 0.15f);
    //static CONSTEXPR float SplitLimit = -1.0f;

    std::vector<vec3> geometricN(m_vertexs.size(), vec3(0));

    std::unordered_multimap<vec3, uint32, std::hash<vec3>> posToVert;

    for (uint32 i = 0; i < m_vertexs.size(); ++i) {
        m_vertexs[i].normal() = vec3(0.0f);
        posToVert.insert(std::make_pair(m_vertexs[i].pos(), i));
    }

    for (TriangleI& t : m_tris) {
        const vec3& p0     = m_vertexs[t.v0].pos();
        const vec3& p1     = m_vertexs[t.v1].pos();
        const vec3& p2     = m_vertexs[t.v2].pos();
        vec3        normal = cross(p1 - p0, (p2 - p0));
        if (normal == vec3(0))
            normal = vec3(0.0f, 1.0f, 0.0f);
        else
            normal = normalize(normal);

        for (int i = 0; i < 3; ++i) {
            vec3& n = geometricN[t.vs[i]];
            if (n == vec3(0)) {
                n = normal;
            } else if (dot(n, normal) < SplitLimit) {
                m_vertexs.push_back(m_vertexs[t.vs[i]]);
                geometricN.push_back(normal);
                t.vs[i] = m_vertexs.size() - 1;
            }
        }
    }

    for (TriangleI& t : m_tris) {
        const vec3& p0     = m_vertexs[t.v0].pos();
        const vec3& p1     = m_vertexs[t.v1].pos();
        const vec3& p2     = m_vertexs[t.v2].pos();
        vec3        normal = cross(p1 - p0, (p2 - p0));
        vec3        nN     = normalize(normal);

        for (int i = 0; i < 3; ++i) {
            auto iters = posToVert.equal_range(m_vertexs[t.vs[i]].pos());

            for (auto t = iters.first; t != iters.second; ++t)
                if (dot(geometricN[t->second], nN) >= SplitLimit)
                    m_vertexs[t->second].normal() += normal;
        }
    }

    for (uint32 i = 0; i < m_vertexs.size(); ++i) {
        if (m_vertexs[i].normal() == vec3(0))
            m_vertexs[i].normal() = geometricN[i];
        else
            m_vertexs[i].normal() = normalize(m_vertexs[i].normal());
    }
}

void TriangleMesh::loadResources(const Json& json, const Scene& scene) {
    MeshIO::LoadMeshFromFile(json["file"].get<std::string>(), m_vertexs, m_tris);
    mat4 transformMatrix = getOptional(json, "transform", getIndentifyTransform());
    Json bsdf_json       = json["bsdf"];
    if (bsdf_json.is_array()) {
        for (const auto& subBsdf : bsdf_json) {
            m_bsdfs.push_back(scene.fetchBSDF(subBsdf));
        }
    } else {
        m_bsdfs.push_back(scene.fetchBSDF(bsdf_json));
    }
    if (transformMatrix != mat4()) {
        mat4        transformNormalMat = getTransformNormalMat(transformMatrix);
        std::string s                  = Mat4ToStr(transformNormalMat);
        for (Vertex& vertex : m_vertexs) {
            vertex.pos()    = transformPoint(transformMatrix, vertex.pos());
            vertex.normal() = transformVector(transformNormalMat, vertex.normal());
        }
    }
}

void TriangleMesh::buildRTC() {
    m_scene    = rtcNewScene(EmbreeUtils::getDevice());
    m_geometry = rtcNewGeometry(EmbreeUtils::getDevice(), RTC_GEOMETRY_TYPE_TRIANGLE);
    // set vertices
    float* vb = (float*)rtcSetNewGeometryBuffer(m_geometry, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * sizeof(float), m_vertexs.size());
    for (size_t i = 0; i < m_vertexs.size(); i++) {
        vb[3 * i]     = m_vertexs[i].pos().x;
        vb[3 * i + 1] = m_vertexs[i].pos().y;
        vb[3 * i + 2] = m_vertexs[i].pos().z;
    }

    // set indices
    unsigned* ib = (unsigned*)rtcSetNewGeometryBuffer(m_geometry, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * sizeof(unsigned), m_tris.size());
    for (size_t i = 0; i < m_tris.size(); i++) {
        ib[i * 3]     = m_tris[i].v0;
        ib[i * 3 + 1] = m_tris[i].v1;
        ib[i * 3 + 2] = m_tris[i].v2;
    }
    for (TriangleI& triangleI : m_tris) triangleI.material = clamp(triangleI.material, 0, m_bsdfs.size() - 1);
    rtcCommitGeometry(m_geometry);
    rtcAttachGeometry(m_scene, m_geometry);
    rtcCommitScene(m_scene);
}

bool TriangleMesh::sameBSDF(BSDF* _bsdf) {
    for (auto bsdf : m_bsdfs)
        if (bsdf.get() == _bsdf)
            return true;
    return false;
}