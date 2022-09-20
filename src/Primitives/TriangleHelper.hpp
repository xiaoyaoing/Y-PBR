#pragma  once

#include "Common/math.hpp"
struct TriangleI{
    union {
        struct { uint32 v0, v1, v2; };
        uint32 vs[3];
    };
    int32 material;
};

struct Vertex{
    vec3 _pos, _normal;
    vec2 _uv;

    const vec3 &normal() const
    {
        return _normal;
    }

    const vec3 &pos() const
    {
        return _pos;
    }

    const vec2 &uv() const
    {
        return _uv;
    }

    vec3 &normal()
    {
        return _normal;
    }

    vec3 &pos()
    {
        return _pos;
    }

    vec2 &uv()
    {
        return _uv;
    }
};




