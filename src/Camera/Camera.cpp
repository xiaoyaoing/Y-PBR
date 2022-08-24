//2022/7/13

#include "Camera.hpp"
#include "../Common/util.hpp"
Camera::Camera(const nlohmann::json & c) {
//    eye = c.at("eye");
//    from_json(c.at("eye"),eye);
//    focal_length = c.at("focal_length").get<Float>() / Float(1000);
//    sensor_width = c.at("sensor_width").get<Float>() / Float(1000);

    _fovDeg=c["fov"];
    _res = c["resolution"];

    if( c.contains("transform")){
        const auto& transform = c["transform"];
        from_json(c["transform"],_transform);
        mat4 & transformMatrix = _transform.matrix;
        _pos    = vec3(transformMatrix[0][3],transformMatrix[1][3],transformMatrix[2][3]);
        _lookAt = _pos+vec3(transformMatrix[0][2],transformMatrix[1][2],transformMatrix[2][2]);
        _up     = vec3(transformMatrix[0][1],transformMatrix[1][1],transformMatrix[2][1]);

        containsAndGet(transform,"lookAt",_lookAt);
        containsAndGet(transform,"up",_up);

        transformMatrix[0][0]=-transformMatrix[0][0];
        transformMatrix[1][0]=-transformMatrix[1][0];
        transformMatrix[2][0]=-transformMatrix[2][0];
//        _ratio = _res.y()/float(_res.x());
//        _pixelSize = Vec2f(1.0f/_res.x(), 1.0f/_res.y());
    }
    spdlog::info("transform Matrix {0}", Mat4ToStr(_transform.matrix));
    preCompute();
//    vec3 look_at=c.at("look_at");
//    lookAt(look_at);
}

void Camera::preCompute( ) {
    _fovRad = Angle::degToRad(_fovDeg);
    _planeDist = 1.0f/std::tan(_fovRad*0.5f);
    _ratio = Float(_res.y)/Float(_res.x);
    _pixelSize = vec2(1.0/_res.x, 1.0/_res.y);
//    float planeArea = (2.0f/_planeDist)*(2.0f*_ratio/_planeDist);
//    _invPlaneArea = 1.0f/planeArea;
}

void Camera::sampleRay(size_t x, size_t y , Ray & ray, vec2 sample) const {

    vec3 localD = normalize(vec3(
            - 1.0 + ( x + 0.5f + sample.x ) * 2.0f * _pixelSize.x,
            _ratio - ( y + 0.5f + sample.y ) * 2.0f * _pixelSize.x,
            _planeDist
    ));

    vec3 d =_transform.matrix * vec4(localD,0);
    d = normalize(d);

    ray = Ray(_pos,d);
}


void Camera::samplePixel(size_t x, size_t y) {

}

void Camera::renderImage(Integrator *integrator, Sampler *sampler) const {

}

//void Camera::lookAt(const vec3& p)
//{
//    forward = normalize(p - eye);
//    left = cross({0.0, 1.0, 0.0}, forward);
//    left = length(left) < Constant::EPSILON ? vec3(-1.0, 0.0, 0.0) : normalize(left);
//    up = normalize(glm::cross(forward, left));
//}








