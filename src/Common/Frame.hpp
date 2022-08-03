#include "math.hpp"
#include "util.hpp"

static void coordinateSystem(const vec3 &a, vec3  &b, vec3  &c) {
    if (std::abs(a.x) > std::abs(a.y)) {
        float invLen = 1.0f / std::sqrt(a.x * a.x + a.z * a.z);
        c = vec3 (a.z * invLen, 0.0f, -a.x * invLen);
    } else {
        float invLen = 1.0f / std::sqrt(a.y * a.y + a.z * a.z);
        c = vec3(0.0f, a.z * invLen, -a.y * invLen);
    }
    vec3  _a(a.x,a.y,a.z);
    b = cross(_a,c);
}

struct Frame {
    vec3 s, t;
    vec3  n;

    /// Default constructor -- performs no initialization!
    Frame() { }

    /// Given a normal and tangent vectors, construct a new coordinate frame
    Frame(const vec3  &s, const vec3 &t, const vec3  &n)
            : s(s), t(t), n(n) { }

    /// Construct a new coordinate frame from a single vector
    Frame(const vec3 n) : n(n){
        coordinateSystem(n, s, t);
    }

    /// Convert from world coordinates to local coordinates
    vec3 toLocal(const vec3 &v) const {
        return vec3(
                dot(v,s), dot(v,t), dot(v,n)
        );
    }

    /// Convert from local coordinates to world coordinates
    vec3 toWorld(const vec3 &v) const {
        return s * v.x + t * v.y   + n * v.z ;
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the cosine of the angle between the normal and v */
    static float cosTheta(const vec3 &v) {
        return v.z;
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the sine of the angle between the normal and v */
    static float sinTheta(const vec3 &v) {
        float temp = sinTheta2(v);
        if (temp <= 0.0f)
            return 0.0f;
        return std::sqrt(temp);
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the tangent of the angle between the normal and v */
    static float tanTheta(const vec3 &v) {
        float temp = 1 - v.z*v.z;
        if (temp <= 0.0f)
            return 0.0f;
        return std::sqrt(temp) / v.z;
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the Squared  sine of the angle between the normal and v */
    static float sinTheta2(const vec3 &v) {
        return 1.0f - v.z * v.z;
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the sine of the phi parameter in spherical coordinates */
    static float sinPhi(const vec3 &v) {
        float sinTheta = Frame::sinTheta(v);
        if (sinTheta == 0.0f)
            return 1.0f;
        return clamp(v.y / sinTheta, -1.0f, 1.0f);
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the cosine of the phi parameter in spherical coordinates */
    static float cosPhi(const vec3 &v) {
        float sinTheta = Frame::sinTheta(v);
        if (sinTheta == 0.0f)
            return 1.0f;
        return clamp(v.x / sinTheta, -1.0f, 1.0f);
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the Squared  sine of the phi parameter in  spherical
     * coordinates */
    static float sinPhi2(const vec3 &v) {
        return clamp(v.y * v.y / sinTheta2(v), 0.0f, 1.0f);
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the Squared  cosine of the phi parameter in  spherical
     * coordinates */
    static float cosPhi2(const vec3 &v) {
        return clamp(v.x * v.x / sinTheta2(v), 0.0f, 1.0f);
    }

    /// Equality test
    bool operator==(const Frame &frame) const {
        return frame.s == s && frame.t == t && frame.n == n;
    }

    /// Inequality test
    bool operator!=(const Frame &frame) const {
        return !operator==(frame);
    }


};