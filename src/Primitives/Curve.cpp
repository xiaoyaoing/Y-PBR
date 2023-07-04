#include "Curve.hpp"
#include "BSpline.hpp"
#include "Sampler/UniformSampler.h"
#include "scene.hpp"
#include <iostream>


static Bounds3 curveBox(const vec4 &q0, const vec4 &q1, const vec4 &q2) {
    vec2 xMinMax(BSpline::quadraticMinMax(q0.x, q1.x, q2.x));
    vec2 yMinMax(BSpline::quadraticMinMax(q0.y, q1.y, q2.y));
    vec2 zMinMax(BSpline::quadraticMinMax(q0.z, q1.z, q2.z));
    float maxW = std::max(q0.w, std::max(q1.w, q2.w));
    return Bounds3(
            vec3(xMinMax.x, yMinMax.x, zMinMax.x) - maxW,
            vec3(xMinMax.y, yMinMax.y, zMinMax.y) + maxW
    );
}

struct CurveIntersection {
    uint32 curveP0;
    float t;
    vec2 uv;
    float w;
};

struct StackNode {
    vec4 p0, p1;
    float tMin, tMax;
    int depth;

    void set(float tMin_, float tMax_, vec4 p0_, vec4 p1_, int depth_) {
        p0 = p0_;
        p1 = p1_;
        tMin = tMin_;
        tMax = tMax_;
        depth = depth_;
    }
};

static vec4 project(const vec3 &o, const vec3 &lx, const vec3 &ly, const vec3 &lz, const vec4 &q) {
    vec3 p(vec3(q) - o);
    return vec4(dot(lx, p), dot(ly, p), dot(lz, p), q.w);
}


static inline void intersectHalfCylinder(StackNode node, float tMin,
                                         float &tMax, CurveIntersection &isect) {

    vec2 v = (node.p1 - node.p0);
    float lengthSq = length2(v);
    float invLengthSq = 1.0f / lengthSq;
    float invLength = std::sqrt(invLengthSq);
    float segmentT = -(dot(vec2(node.p0.x, node.p0.y), v)) * invLengthSq;
    float signedUnnormalized = node.p0.x * v.y - node.p0.y * v.x;
    float distance = std::fabs(signedUnnormalized) * invLength;

    float width = node.p0.w * (1.0f - segmentT) + node.p1.w * segmentT;
    if (distance > width)
        return;

    float depth = node.p0.z * (1.0f - segmentT) + node.p1.z * segmentT;
    float dz = node.p1.z - node.p0.z;
    float ySq = sqr(width) - sqr(distance);
    float lSq = ySq * (1.0f + dz * dz * invLengthSq);
    float deltaT = std::sqrt(std::max(lSq, 0.0f));
    float t0 = depth - deltaT;

    vec3 v3(node.p0 - node.p1);
    lengthSq = length2(v3);
    segmentT = dot((vec3(node.p0.x, node.p0.y, node.p0.z - t0)), v3) / lengthSq;
    if (segmentT < 0.0f || t0 >= tMax || t0 <= tMin) {
        // Enable this for two-sided cylinders (for transmissive BSDFs)
        // Not really recommended (self-intersections), so disabled for now

//        t0 = depth + deltaT;
//        segmentT -= 2.0f*deltaT*v3.z/lengthSq;
//        if (segmentT < 0.0f || t0 >= closestDepth || t0 <= tMin)
        return;
    }

    float newT = segmentT * (node.tMax - node.tMin) + node.tMin;

    if (newT >= 0.0f && newT <= 1.0f) {
        isect.uv = vec2(newT, 0.5f + 0.5f * distance / width);
        isect.t = t0;
        isect.w = width;
        tMax = t0;
    }
}

template<typename T>
static inline void precomputeBSplineCoefficients(T &p0, T &p1, T &p2) {
    T q0 = (0.5f * p0 - p1 + 0.5f * p2);
    T q1 = (p1 - p0);
    T q2 = 0.5f * (p0 + p1);
    p0 = q0;
    p1 = q1;
    p2 = q2;
}

static bool pointOnSpline(vec4 q0, vec4 q1, vec4 q2,
                          float tMin, float tMax, CurveIntersection &isect,
                          vec3 n0, vec3 n1, vec3 n2) {
    const int MaxDepth = 5;

    StackNode stackBuf[MaxDepth];
    StackNode *stack = &stackBuf[0];

    precomputeBSplineCoefficients(q0, q1, q2);

    vec4 tFlat = -q1 * 0.5f / q0;
    vec2 xyFlat = (q0 * tFlat * tFlat + q1 * tFlat + q2);
    float xFlat = xyFlat.x, yFlat = xyFlat.y;

    StackNode cur{
            q2,
            q0 + q1 + q2,
            0.0f, 1.0f, 0
    };
    float closestDepth = tMax;

    while (true) {
        vec2 pMin = min(cur.p0, cur.p1);
        vec2 pMax = max(cur.p0, cur.p1);
        if (tFlat.x > cur.tMin && tFlat.x < cur.tMax) {
            pMin.x = std::min(pMin.x, xFlat);
            pMax.x = std::max(pMax.x, xFlat);
        }
        if (tFlat.y > cur.tMin && tFlat.y < cur.tMax) {
            pMin.y = std::min(pMin.y, yFlat);
            pMax.y = std::max(pMax.y, yFlat);
        }

        float maxWidth = std::max(cur.p0.w, cur.p1.w);
        if (pMin.x <= maxWidth && pMin.y <= maxWidth && pMax.x >= -maxWidth && pMax.y >= -maxWidth) {
            if (cur.depth >= MaxDepth) {
                intersectHalfCylinder(cur, tMin, closestDepth, isect);
            } else {
                float splitT = (cur.tMin + cur.tMax) * 0.5f;
                vec4 qSplit = q0 * (splitT * splitT) + q1 * splitT + q2;

                if (cur.p0.z < qSplit.z) {
                    (*stack++).set(splitT, cur.tMax, qSplit, cur.p1, cur.depth + 1);
                    cur.set(cur.tMin, splitT, cur.p0, qSplit, cur.depth + 1);
                } else {
                    (*stack++).set(cur.tMin, splitT, cur.p0, qSplit, cur.depth + 1);
                    cur.set(splitT, cur.tMax, qSplit, cur.p1, cur.depth + 1);
                }
                continue;
            }
        }
        do {
            if (stack == stackBuf)
                return closestDepth < tMax;
            cur = *--stack;
        } while (std::min(cur.p0.z - cur.p0.w, cur.p1.z - cur.p1.w) > closestDepth);
    }

    return false;
}

void Curve::computeArea() {
    Primitive::computeArea();
}

void Curve::computeBoundingBox() {
    BB_ = Bounds3();
    for (size_t i = 2; i < _nodeData.size(); ++i)
        BB_ = Union(BB_, curveBox(_nodeData[i - 2], _nodeData[i - 1], _nodeData[i]));
    // BB_  = Bounds3(vec3(-1000),vec3(1000));
}

void Curve::transform(const mat4 &T) {
}

static bool useEmbree = true;

inline vec3 derivBezier(const std::vector<vec4> &positions, const unsigned int primID, const float t) {
    vec3 p00, p01, p02, p03;
    p00 = positions[primID];
    p01 = positions[primID + 1];
    p02 = positions[primID + 2];
    p03 = positions[primID + 3];
    const float t0 = 1.0f - t, t1 = t;
    const vec3 p10 = p00 * t0 + p01 * t1;
    const vec3 p11 = p01 * t0 + p02 * t1;
    const vec3 p12 = p02 * t0 + p03 * t1;
    const vec3 p20 = p10 * t0 + p11 * t1;
    const vec3 p21 = p11 * t0 + p12 * t1;
    //return p20 * t0 + p21 * t1;
    return vec3(3.0f * (p21 - p20));
}

inline vec3 derivSpline(const std::vector<vec4> &positions, const unsigned int primID, const float t) {
    vec3 p00, p01, p02, p03;
    p00 = positions[primID];
    p01 = positions[primID + 1];
    p02 = positions[primID + 2];
    p03 = positions[primID + 3];
    //return p20 * t0 + p21 * t1;
    const float t0 = 1.0f - t, t1 = t;
    const float n0 = -0.5f * t1 * t1;
    const float n1 = -0.5f * t0 * t0 - 2.0f * (t0 * t1);
    const float n2 = 0.5f * t1 * t1 + 2.0f * (t1 * t0);
    const float n3 = 0.5f * t0 * t0;
    return vec3(n0 * p00 + n1 * p01 + n2 * p02 + n3 * p03);
}


std::optional<Intersection> Curve::intersect(Ray &ray) const {
    // return std::nullopt;
    if (useEmbree) {
        RTCRayHit rayHit;
        EmbreeUtils::convertRay(&ray, &rayHit);

        RTCIntersectContext context;
        rtcInitIntersectContext(&context);

        rtcIntersect1(m_scene, &context, &rayHit);

        if (rayHit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {

            float P[4], P1[4], P2[4];
            rtcInterpolate1(rtcGetGeometry(m_scene, rayHit.hit.geomID), rayHit.hit.primID, rayHit.hit.u, rayHit.hit.v,
                            RTC_BUFFER_TYPE_VERTEX, 0, P, P1, P2, 4);

            ray.farT = rayHit.ray.tfar;
            Intersection its;
            its.p = ray.operator()(ray.farT);
            its.bsdf = bsdf.get();
            its.primitive = this;
            its.w = -ray.d;
            uint32 p0 = _indices[rayHit.hit.primID];
            float t = rayHit.hit.u;
          //  p0 = 5;t=0.5;
            vec3 tangent = normalize(
                    BSpline::quadraticDeriv(_nodeData[p0], _nodeData[p0 + 1], _nodeData[p0 + 2],  t));

            vec3 point = BSpline::quadratic(_nodeData[p0], _nodeData[p0 + 1], _nodeData[p0 + 2], _nodeData[p0 + 3], t);
            its.Ng = its.Ns = normalize((its.w - tangent * dot(tangent, its.w)));
            its.tangent = {vec3(tangent)};
            its.uv = {rayHit.hit.u, rayHit.hit.v};
//            its.Ng = its.tangent.value();
//            its.Ng = vec3(rayHit.hit.Ng_x,rayHit.hit.Ng_y,rayHit.hit.Ng_z);
//            its.Ng = vec3(float(p0)/_nodeData.size());
            return {its};
//            its.Ng = its.Ns = normalize(v);
//            its.Ng = vec3(t)*2.f - 1.f;
//            return {its};
//
//            vec3 dp = derivSpline(_nodeData,p0,t);
//            auto Tx = normalize(dp);
//            auto Ty = normalize(cross(-(ray.d),Tx));
//
//            its.Ng = its.Ns = normalize(cross(Ty,Tx));
//            its.Ng = its.Ns = vec3(rayHit.hit.Ng_x,rayHit.hit.Ng_y,rayHit.hit.Ng_z);
//            return {its};
        }
        return std::nullopt;
    }

    CurveIntersection curveIts;
    bool didIntersect = false;
    auto len = length(ray.d);
    vec3 o(ray.o), lz(ray.d);
    float d = std::sqrt(lz.x * lz.x + lz.z * lz.z);
    vec3 lx, ly;
    if (d == 0.0f) {
        lx = vec3(1.0f, 0.0f, 0.0f);
        ly = vec3(0.0f, 0.0f, -lz.y);
    } else {
        lx = vec3(lz.z / d, 0.0f, -lz.x / d);
        ly = vec3(lx.z * lz.y, d, -lz.y * lx.x);
    }


    _bvh->trace(ray, [&](Ray &ray, uint32 id) {

        vec4 q0(project(o, lx, ly, lz, _nodeData[id - 2]));
        vec4 q1(project(o, lx, ly, lz, _nodeData[id - 1]));
        vec4 q2(project(o, lx, ly, lz, _nodeData[id - 0]));
        auto t = curveBox(_nodeData[id - 2], _nodeData[id - 1], _nodeData[id - 0]);
        vec3 invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
        int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
        bool temp = t.IntersectP(ray, invDir, dirIsNeg);
        vec3 n0, n1, n2;

        if (pointOnSpline(q0, q1, q2, ray.nearT, ray.farT, curveIts, n0, n1, n2)) {
            ray.farT = curveIts.t;
            curveIts.curveP0 = id - 2;
            didIntersect = true;
        }
    });

    if (!didIntersect) {
        return std::nullopt;
    }

    Intersection its;
    its.w = -ray.d;
    uint32 p0 = curveIts.curveP0;
    float t = curveIts.uv.x;

    vec3 tangent = normalize(BSpline::quadraticDeriv(_nodeData[p0], _nodeData[p0 + 1], _nodeData[p0 + 2], t));

    if (_mode == MODE_RIBBON) {
        vec3 normal = BSpline::quadratic(_nodeNormals[p0], _nodeNormals[p0 + 1], _nodeNormals[p0 + 2], t);
        its.Ng = its.Ns = normalize(dot(tangent * tangent, normal) - normal);
    } else if (_mode == MODE_BCSDF_CYLINDER) {
        its.Ng = its.Ns = normalize(-its.w - dot(tangent * tangent, -its.w));
    } else if (_mode == MODE_HALF_CYLINDER || _mode == MODE_CYLINDER) {
        vec3 point = BSpline::quadratic(_nodeData[p0], _nodeData[p0 + 1], _nodeData[p0 + 2], t);
        vec3 localP = its.p - point;
        localP -= tangent * (dot(localP, tangent));
        its.Ng = its.Ns = normalize(localP);
    }
    its.uv = curveIts.uv;
    its.primitive = this;
    its.bsdf = bsdf.get();

    return {its};

//    if (_mode == MODE_CYLINDER)
//        its.epsilon = max(its.epsilon, 0.1f*isect.w);
//    else
//        its.epsilon = max(its.epsilon, 0.01f*isect.w);
}

bool Curve::occluded(const Ray &ray) const {
    Ray _ray(ray);
    return intersect(_ray).has_value();
}

Curve::Curve(const Json &json, Scene &scene) : Primitive(json) {
    CurveIO::CurveData data{&_curveEnds, &_nodeData, &_nodeColor, &_nodeNormals};
    CurveIO::load(json["file"], data);
    _curveCount = _curveEnds.size();

    mat4 toWorld = getOptional(json, "transform", getIndentifyTransform());

    _overrideThickness = containsAndGet(json, "curve_thickness", _curveThickness);
    bool tapper = getOptional(json, "curve_taper", false);
    if (_overrideThickness || tapper) {
        for (uint32 i = 0; i < _curveCount; ++i) {
            uint32 start = i ? _curveEnds[i - 1] : 0;
            for (uint32 t = start; t < _curveEnds[i]; ++t) {
                float thickness = _overrideThickness ? _curveThickness : _nodeData[t].w;
                if (tapper)
                    thickness *= 1.0f - (t - start - 0.5f) / (_curveEnds[i] - start - 1);
                _nodeData[t].w = thickness;
            }
        }
    }

    vec3 widthScaleVec = extractScaleVec(toWorld);
    Float widthScale = (widthScaleVec.x + widthScaleVec.y + widthScaleVec.z) / 3;
    for (vec4 &nodeData: _nodeData) {
        vec3 newP = transformPoint(toWorld, vec3(nodeData.x, nodeData.y, nodeData.z));
        nodeData = vec4(newP.x, newP.y, newP.z, nodeData.w * widthScale);
    }
    _subSample = 0;

    std::vector<std::shared_ptr<Primitive>> prims;
    UniformSampler rand(_curveCount);
    for (uint32 i = 0; i < _curveCount; i++) {
        uint32 start = 0;
        if (i > 0) {
            start = _curveEnds[i - 1];
        }
        for (uint32 t = start + 2; t < _curveEnds[i]; ++t) {
            const vec4 &p0 = _nodeData[t - 2];
            const vec4 &p1 = _nodeData[t - 1];
            const vec4 &p2 = _nodeData[t - 0];

            if (_subSample > 0.0f && rand.getNext1D() < _subSample)
                continue;
//            std::shared_ptr<BSDF> bsdf = scene.fetchBSDF(getOptional(json, "bsdf", Json()));
//            //   prim->setBSDF(bsdf);
//            //  scene.AddPrimitive(prim);
//            );
        }

    }


    std::cout << "Building Hair BVh Start" << std::endl;
   // if (!useEmbree) _bvh.reset(new BVHAccel(prims, BVHAccel::SplitMethod::PBRTSAH, 2));
    std::cout << "Building Hair BVh End" << std::endl;
    computeBoundingBox();

    m_scene = rtcNewScene(EmbreeUtils::getDevice());
    m_geometry = rtcNewGeometry(EmbreeUtils::getDevice(), RTC_GEOMETRY_TYPE_FLAT_BSPLINE_CURVE);
    float *vb = (float *) rtcSetNewGeometryBuffer(m_geometry, RTC_BUFFER_TYPE_VERTEX, 0,
                                                  RTC_FORMAT_FLOAT4,
                                                  4 * sizeof(float), _nodeData.size());
    for (size_t i = 0; i < _nodeData.size(); i++) {
        vb[4 * i] = _nodeData[i].x;
        vb[4 * i + 1] = _nodeData[i].y;
        vb[4 * i + 2] = _nodeData[i].z;
        vb[4 * i + 3] = _nodeData[i].w;
    }


//    rtcSetSharedGeometryBuffer
//    (m_geometry,RTC_BUFFER_TYPE_VERTEX,0,RTC_FORMAT_FLOAT4,&(_nodeData[0]),sizeof(vec4),0,_nodeData.size());
//    rtcSetSharedGeometryBuffer
//    (m_geometry,RTC_BUFFER_TYPE_INDEX,0,RTC_FORMAT_UINT,&(_curveEnds[0]),sizeof(uint32),0,_curveEnds.size());
    uint32 start = 0;
    for (int i = 0; i < _curveCount; i++) {
        if (i > 0) {
            start = _curveEnds[i - 1];
        }
        for (uint32 t = start + 3; t < _curveEnds[i]; ++t) {
            _indices.push_back(t - 3);
        }
    }

    unsigned *ib = (unsigned *) rtcSetNewGeometryBuffer(m_geometry, RTC_BUFFER_TYPE_INDEX, 0,
                                                        RTC_FORMAT_UINT,
                                                        sizeof(unsigned), _indices.size());
//    for ( size_t i = 0 ; i < _curveEnds.size() ; i ++ ) {
//        ib[i] = _curveEnds[i];
//    }
    memcpy(ib, _indices.data(), _indices.size());
    rtcCommitGeometry(m_geometry);
    rtcAttachGeometry(m_scene, m_geometry);
    rtcCommitScene(m_scene);
}

Frame Curve::setTangentFrame(const Intersection *its) const {
    vec3 T, B, N;
    N = its->Ns;
    B = its->tangent.value();
    T = cross(B, N);
    T = normalize(T - dot(T, N) * N);
    B = cross(N, T);

    if (length(T) < Constant::EPSILON || length(B) < Constant::EPSILON) {

    }

    return Frame(T, B, N);
}

//std::optional<Intersection> CurveI::intersect(Ray &ray) const {
//    CurveIntersection curveIts;
//    bool didIntersect = false;
//    vec3 o(ray.o), lz(ray.d);
//    float d = std::sqrt(lz.x * lz.x + lz.z * lz.z);
//    vec3 lx, ly;
//    if (d == 0.0f) {
//        lx = vec3(1.0f, 0.0f, 0.0f);
//        ly = vec3(0.0f, 0.0f, -lz.y);
//    } else {
//        lx = vec3(lz.z / d, 0.0f, -lz.x / d);
//        ly = vec3(lx.z * lz.y, d, -lz.y * lx.x);
//    }
//
//    vec4 q0(project(o, lx, ly, lz, (*_nodeData)[id - 2]));
//    vec4 q1(project(o, lx, ly, lz, (*_nodeData)[id - 1]));
//    vec4 q2(project(o, lx, ly, lz, (*_nodeData)[id - 0]));
//    vec3 invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
//    int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
//    vec3 n0, n1, n2;
//
//    if (pointOnSpline(q0, q1, q2, ray.nearT, ray.farT, curveIts, n0, n1, n2)) {
//        ray.farT = curveIts.t;
//        curveIts.curveP0 = id - 2;
//        didIntersect = true;
//    }
//
//    if (!didIntersect) {
//        return std::nullopt;
//    }
//
//    Intersection its;
//    its.w = -ray.d;
//    its.p = ray.o + ray.d * ray.farT;
//    uint32 p0 = curveIts.curveP0;
//    float t = curveIts.uv.x;
//    t = 0.5;
//    vec3 tangent = normalize(BSpline::quadraticDeriv((*_nodeData)[p0], (*_nodeData)[p0 + 1], (*_nodeData)[p0 + 2], t));
//    {
//        vec3 point = BSpline::quadratic((*_nodeData)[p0], (*_nodeData)[p0 + 1], (*_nodeData)[p0 + 2], t);
//        vec3 localP = its.p - point;
//        localP -= tangent * (dot(localP, tangent));
//        its.Ng = its.Ns = normalize(localP);
//    }
//    its.uv = curveIts.uv;
//    its.primitive = this;
//    its.bsdf = bsdf.get();
//
//    return {its};
//}
//
//bool CurveI::occluded(const Ray &ray) const {
//    Ray _ray(ray);
//    return intersect(_ray).has_value();
//}
//
//void CurveI::computeBoundingBox() {
//    BB_ = curveBox((*_nodeData)[id - 2], (*_nodeData)[id - 1], (*_nodeData)[id - 0]);
//}
