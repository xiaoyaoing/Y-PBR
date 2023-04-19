//2022/7/13

#include "Camera.hpp"
#include "../Common/Json.hpp"

Camera::Camera(const Json &c) {

    _fovDeg = getOptional(c, "fov", 45);
    _res = getOptional(c, "resolution", ivec2(512, 512));
    image = std::make_unique<Image>(_res, getOptional(c, "output", std::string("image")));
    if (c.contains("transform")) {
        const auto &transformJson = c["transform"];
        _cameraToWorld = c["transform"];
        _pos = vec3(_cameraToWorld[0][3], _cameraToWorld[1][3], _cameraToWorld[2][3]);
        _lookAt = _pos + vec3(_cameraToWorld[0][2], _cameraToWorld[1][2], _cameraToWorld[2][2]);
        _up = vec3(_cameraToWorld[0][1], _cameraToWorld[1][1], _cameraToWorld[2][1]);
        containsAndGet(transformJson, "lookAt", _lookAt);
        containsAndGet(transformJson, "up", _up);
        _cameraToWorld[0][0] = -_cameraToWorld[0][0];
        _cameraToWorld[1][0] = -_cameraToWorld[1][0];
        _cameraToWorld[2][0] = -_cameraToWorld[2][0];
    }
    _toLocal = glm::inverse(_cameraToWorld);

    preCompute();
}

void Camera::preCompute() {
    _fovRad = Angle::degToRad(_fovDeg);
    _planeDist = 1.0f / std::tan(_fovRad * 0.5f);
    _ratio = Float(_res.y) / Float(_res.x);
    _pixelSize = vec2(1.0 / _res.x, 1.0 / _res.y);

    _rasterToCamera = getTransFormMatrix(vec3(-1, _ratio, _planeDist), vec3(2.f / _res.x, -2.f / _res.x, 0), vec3());
    _cameraToRaster = glm::inverse(_rasterToCamera);
}

Ray Camera::sampleRay(size_t x, size_t y, vec2 sample) const {
    auto localD = normalize(transformPoint(_rasterToCamera, vec3(x + sample.x, y + sample.y, 0)));
//    vec3 localD = normalize(vec3(
//            - 1.0 + ( x + 0.5f + sample.x ) * 2.0f * _pixelSize.x,
//            _ratio - ( y + 0.5f + sample.y ) * 2.0f * _pixelSize.x,
//            _planeDist
//    ));
    vec3 d = transformVector(_cameraToWorld, localD);
    d = normalize(d);
    return Ray(_pos, d);
}


vec2 Camera::inverse(const vec3 &pos) const {
    vec3 worldDir = normalize(pos - _pos);
    vec3 localD = transformVector(_toLocal, worldDir);
    localD *= _planeDist / localD.z;
    return {(localD.x + 1) / (2.f * _pixelSize.x), (-localD.y + _ratio) / (2 * _pixelSize.x)};
}

static Float factor(Float a, Float b, Float c) {
    if (a < 0 && b < 0) {
        a = -a;
        b = -b;
    }
    return (c - a) / (b - a);
}

void Camera::drawLine(int x0, int x1, int y0, int y1, const Spectrum &color) {
    int dy = y1 - y0;
    int dx = x1 - x0;
    int dy2 = 2 * dy;
    int dx2 = 2 * dx;

    if (dy == 0) {
        if (x0 > x1) std::swap(x0, x1);
        for (int x = x0; x < x1; x++)
            image->addPixel(x, y0, color);
        return;
    }
    if (dx == 0) {
        if (y0 > y1) std::swap(y0, y1);
        for (int y = y0; y < y1; y++)
            image->addPixel(x0, y, color);
        return;
    }

    Float m = abs(dy / dx);
    if (dx > 0 && dy > 0 && m < 1) {
        Float pi = dy2 - dx;
        for (int x = x0, y = y0; x < x1; x++) {
            if (pi < 0)
                pi += dy2;
            else {
                pi += dy2 - dx2;
                y += 1;
            }
            image->addPixel(x, y, color * factor(x0, x1, x));
        }
    } else if (dx > 0 && dy > 0 && m >= 1) {
        Float pi = dx2 - dy;
        for (int x = x0, y = y0; y < y1; y++) {
            if (pi < 0)
                pi += dx2;
            else {
                pi += dx2 - dy2;
                x += 1;
            }
            image->addPixel(x, y, color * factor(y0, y1, y));
        }
    } else if (dx > 0 && dy < 0 && m < 1) {
        dy2 = -dy2;
        dy = -dy;
        Float pi = dy2 - dx;
        for (int x = x0, y = -y0; x < x1; x++) {
            if (pi < 0)
                pi += dy2;
            else {
                pi += dy2 - dx2;
                y += 1;
            }
            image->addPixel(x, -y, color * factor(x0, x1, x));
        }
    } else if (dx > 0 && dy < 0 && m >= 1) {
        dy2 = -dy2;
        dy = -dy;
        Float pi = dx - dy2;
        for (int x = x0, y = -y0; y < -y1; y++) {
            if (pi < 0)
                pi += dx2;
            else {
                pi += dx2 - dy2;
                x += 1;
            }
            image->addPixel(x, -y, color * factor(y0, y1, y));
        }
    } else if (dx < 0 && dy > 0 && m < 1) {
        dx2 = -dx2;
        dx = dx;
        Float pi = dy2 - dx;
        for (int x = -x0, y = y0; x < -x1; x++) {
            if (pi < 0)
                pi += dy2;
            else {
                pi += dy2 - dx2;
                y += 1;
            }
            image->addPixel(-x, y, color * factor(x0, x1, x));
        }
    } else if (dx < 0 && dy > 0 && m >= 1) {
        dx2 = -dx2;
        dx = dx;
        Float pi = dx - dy2;
        for (int x = x0, y = -y0; y < -y1; y++) {
            if (pi < 0)
                pi += dx2;
            else {
                pi += dx2 - dy2;
                x += 1;
            }
            image->addPixel(-x, y, color * factor(y0, y1, y));
        }
    } else if (dx < 0 && dy < 0 && m < 1) {
        dx2 = -dx2;
        dx = -dx;
        dy2 = -dy2;
        dy = -dy;
        Float pi = dy2 - dx;
        for (int x = -x0, y = -y0; x < -x1; x++) {
            if (pi < 0)
                pi += dy2;
            else {
                pi += dy2 - dx2;
                y += 1;
            }
            image->addPixel(-x, -y, color * factor(x0, x1, x));
        }
    } else if (dx < 0 && dy < 0 && m >= 1) {
        dx2 = -dx2;
        dx = -dx;
        dy2 = -dy2;
        dy = -dy;
        Float pi = dx - dy2;
        for (int x = -x0, y = -y0; y < -y1; y++) {
            if (pi < 0)
                pi += dx2;
            else {
                pi += dx2 - dy2;
                x += 1;
            }
            image->addPixel(-x, -y, color * factor(y0, y1, y));
        }
    }
}

static bool draw = false;
static int count = 0;

void Camera::drawLine(const vec3 &begin, const vec3 &end, const Spectrum &color) {
    ivec2 screenBegin = inverse(begin);
    ivec2 screenEnd = inverse(end);
    screenBegin = clamp(screenBegin, ivec2(0), image->resoulation());
    screenEnd = clamp(screenEnd, ivec2(0), image->resoulation());
    int x0 = screenBegin.x, x1 = screenEnd.x;
    int y0 = screenBegin.y, y1 = screenEnd.y;
    drawLine(x0, x1, y0, y1, color);
    draw = true;
}

Spectrum Camera::evalDirection(vec3 p, ivec2 *pRaster, Float *pdf) const {
    auto cameraP = transformPoint(_cameraToWorld, vec3(0, 0, 0));
    auto ray = Ray(cameraP, direction(cameraP, p));
    auto sw = rayWeight(ray, pRaster);
    if (!isBlack(sw))
        *pdf = distance2(p, cameraP) / dot(ray.d, transformVector(_cameraToWorld, vec3(0, 0, 1)));
    return sw;
}

//默认是从相机那一点出发
Spectrum Camera::rayWeight(const Ray &ray, ivec2 *pRaster) const {
    auto localO = transformPoint(_toLocal, ray.o);
    auto localD = transformVector(_toLocal, ray.d);
    auto cosTheta = dot(localD, vec3(0, 0, 1));
    auto cos2Theta = cosTheta * cosTheta;
    if (cosTheta < 0) {
        return Spectrum(0);
    }
    auto p = ray(_planeDist / cosTheta);
    auto pRaster3d = transformVector(_cameraToRaster, transformPoint(_toLocal, p));
    if (pRaster3d.x < _res.x && pRaster3d.y < _res.y) {
        if (pRaster)
            *pRaster = ivec2(pRaster3d.x, pRaster3d.y);
        return Spectrum(1.f / (A * cos2Theta * cos2Theta));
    } else {
        return Spectrum(0);
    }
}

//pinole sample
bool Camera::samplePosition(ivec2 point, vec2 sample, PositionSample &pSample) const {
    pSample.p = _pos;
    pSample.weight = Spectrum(1);
    pSample.pdf = 1;
    pSample.normal = transformVector(_cameraToWorld, vec3(0, 0, 1));
    return true;
}

PositionAndDirectionSample Camera::sampleRay(ivec2 point, vec2 /*posSample*/, vec2 dirSample) const {
    PositionAndDirectionSample result;
    auto localD = normalize(transformPoint(_rasterToCamera, vec3(point.x + dirSample.x, point.y + dirSample.y, 0)));
    result.posPdf = 1;
    result.dirPdf = _invPlaneArea / (localD.z * localD.z * localD.z);
    result.weight = Spectrum(1);
    result.n = transformVector(_cameraToWorld, vec3(0, 0, 1));
    auto worldD = transformVector(_cameraToWorld,localD);
    result.ray = Ray(_pos,worldD);
    return result;
}

PositionAndDirectionSample Camera::sampleLi(vec3 p, ivec2 *pRaster, vec2 sample) const {
    PositionAndDirectionSample result;
    auto cameraP = transformPoint(_cameraToWorld, vec3(0, 0, 0));
    result.ray = Ray(cameraP, direction(cameraP, p));
    result.weight = rayWeight(result.ray, pRaster);
    if (!isBlack(result.weight))
        result.dirPdf = distance2(p, cameraP) / dot(result.ray.d, transformVector(_cameraToWorld, vec3(0, 0, 1)));
    result.n= transformVector(_cameraToWorld, vec3(0, 0, 1));
    result.posPdf = 1;
    result.weight/=(result.dirPdf * result.posPdf);
    return result;
    }








