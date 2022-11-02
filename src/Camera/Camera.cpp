//2022/7/13

#include "Camera.hpp"
#include "../Common/Json.hpp"

Camera::Camera(const Json & c) {

    _fovDeg = getOptional(c, "fov", 45);
    _res = getOptional(c, "resolution", ivec2(512, 512));
    image = std::make_unique < Image >(_res, getOptional(c, "output", std::string("image")));
    if ( c.contains("transform") ) {
        const auto & transformJson = c["transform"];
        _toWorld = c["transform"];
        _pos = vec3(_toWorld[0][3], _toWorld[1][3], _toWorld[2][3]);
        _lookAt = _pos + vec3(_toWorld[0][2], _toWorld[1][2], _toWorld[2][2]);
        _up = vec3(_toWorld[0][1], _toWorld[1][1], _toWorld[2][1]);
        containsAndGet(transformJson, "lookAt", _lookAt);
        containsAndGet(transformJson, "up", _up);
        _toWorld[0][0] = - _toWorld[0][0];
        _toWorld[1][0] = - _toWorld[1][0];
        _toWorld[2][0] = - _toWorld[2][0];
    }
    _toLocal = glm::inverse(_toWorld);
    preCompute();
}

void Camera::preCompute( ) {
    _fovRad = Angle::degToRad(_fovDeg);
    _planeDist = 1.0f / std::tan(_fovRad * 0.5f);
    _ratio = Float(_res.y) / Float(_res.x);
    _pixelSize = vec2(1.0 / _res.x, 1.0 / _res.y);

}

Ray Camera::sampleRay(size_t x, size_t y, vec2 sample) const {
    vec3 localD = normalize(vec3(
            - 1.0 + ( x + 0.5f + sample.x ) * 2.0f * _pixelSize.x,
            _ratio - ( y + 0.5f + sample.y ) * 2.0f * _pixelSize.x,
            _planeDist
    ));
    vec3 d = transformVector(_toWorld, localD);
    d = normalize(d);
    return Ray(_pos, d);
}


vec2 Camera::inverse(const vec3 & pos) const {
    vec3 worldDir = normalize(pos-_pos);
    vec3 localD = transformVector(_toLocal,worldDir) ;
    localD *= _planeDist / localD.z;
    return {(localD.x+1)/(2.f * _pixelSize.x),(-localD.y+_ratio)/(2 * _pixelSize.x)};
}

static Float factor(Float a,Float b,Float c){
    if(a<0 && b<0) {a=-a;b=-b;}
    return (c-a)/(b-a);
}

void Camera::drawLine(int x0,int x1,int y0,int y1, const Spectrum & color) {
    int dy = y1-y0;
    int dx = x1-x0;
    int dy2 = 2*dy;
    int dx2 = 2*dx;

    if(dy==0){
        if(x0>x1) std::swap(x0,x1);
        for(int x = x0;x<x1;x++)
            image->addPixel(x,y0,color);
        return ;
    }
    if(dx==0){
        if(y0>y1) std::swap(y0,y1);
        for(int y = y0;y<y1;y++)
            image->addPixel(x0,y,color);
        return ;
    }

    Float m = abs(dy/dx);
    if(dx> 0 && dy>0 && m<1 ){
        Float pi = dy2 -dx;
        for(int x = x0,y=y0;x<x1;x++){
            if(pi<0)
                pi += dy2;
            else{
                pi += dy2 -dx2;
                y+=1;
            }
            image->addPixel(x,y,color * factor(x0,x1,x)) ;
        }
    }
    else if(dx>0 && dy>0 && m>=1){
        Float pi = dx2 - dy;
        for(int x = x0,y=y0;y<y1;y++){
            if(pi<0)
                pi += dx2;
            else{
                pi += dx2 -dy2;
                x+=1;
            }
            image->addPixel(x,y,color * factor(y0,y1,y));
        }
    }
    else if(dx>0 && dy<0 && m<1){
        dy2 = -dy2; dy = -dy;
        Float pi = dy2 - dx;
        for(int x = x0,y=-y0;x<x1;x++){
            if(pi<0)
                pi += dy2;
            else{
                pi += dy2 -dx2;
                y+=1;
            }
            image->addPixel(x,-y,color * factor(x0,x1,x));
        }
    }else if(dx>0 && dy<0 && m>=1){
        dy2 = -dy2; dy = -dy;
        Float pi = dx -dy2;
        for(int x = x0,y=-y0;y<-y1;y++){
            if(pi<0)
                pi += dx2;
            else{
                pi += dx2 -dy2;
                x+=1;
            }
            image->addPixel(x,-y,color * factor(y0,y1,y) );
        }
    }else if(dx<0 && dy>0 && m<1){
        dx2 = -dx2; dx = dx;
        Float pi = dy2 -dx;
        for(int x = -x0,y=y0;x<-x1;x++){
            if(pi<0)
                pi += dy2;
            else{
                pi += dy2 -dx2;
                y+=1;
            }
            image->addPixel(-x,y,color * factor(x0,x1,x));
        }
    }else if(dx<0 && dy>0 && m>=1){
        dx2 = -dx2; dx = dx;
        Float pi = dx -dy2;
        for(int x = x0,y=-y0;y<-y1;y++){
            if(pi<0)
                pi += dx2;
            else{
                pi += dx2 -dy2;
                x+=1;
            }
            image->addPixel(-x,y,color * factor(y0,y1,y));
        }
    }else if(dx<0 && dy<0 && m<1){
        dx2 = -dx2; dx = -dx;
        dy2 = -dy2; dy = -dy;
        Float pi = dy2 - dx;
        for(int x = -x0,y=-y0;x<-x1;x++){
            if(pi<0)
                pi += dy2;
            else{
                pi += dy2 -dx2;
                y+=1;
            }
            image->addPixel(-x,-y,color * factor(x0,x1,x));
        }
    }else if(dx<0 && dy<0 && m>=1){
        dx2 = -dx2; dx = -dx;
        dy2 = -dy2; dy = -dy;
        Float pi = dx -dy2;
        for(int x = -x0,y= -y0;y<-y1;y++){
            if(pi<0)
                pi += dx2;
            else{
                pi += dx2 -dy2;
                x+=1;
            }
            image->addPixel(-x,-y,color * factor(y0,y1,y));
        }
    }
}

static bool draw = false;
static int count  = 0;
void Camera::drawLine(const vec3 & begin, const vec3 & end, const Spectrum & color) {
    ivec2 screenBegin = inverse(begin);
    ivec2 screenEnd = inverse(end);
    screenBegin = clamp(screenBegin,ivec2(0),image->resoulation());
    screenEnd = clamp(screenEnd,ivec2(0),image->resoulation());
    int x0=screenBegin.x,x1 = screenEnd.x;
    int y0=screenBegin.y,y1 = screenEnd.y;
   // if(co) return ;
    //inverse(end);
    drawLine(x0,x1,y0,y1,color);
//    int xx[] = {0,255};
//    int yy[] = {0,128};
//    for(int i =0 ;i<8;i++){
//        x0 = xx[i &1],x1 = xx[1-i&1];
//        y0 = yy[ i>>1 & 1 ],y1 = yy[1-(i>>1 &1)];
//        if(i>>4 & 1){
//            std::swap(x0,y0);
//            std::swap(x1,y1);
//        }
//        //drawLine(0,128,0,255,color);
//        drawLine(0,128,255,0,color);
//       // drawLine(x0,x1,y0,y1,color);
//    }

    draw = true;
}


//void Camera::lookAt(const vec3& p)
//{
//    forward = normalize(p - eye);
//    left = cross({0.0, 1.0, 0.0}, forward);
//    left = length(left) < Constant::EPSILON ? vec3(-1.0, 0.0, 0.0) : normalize(left);
//    up = normalize(glm::cross(forward, left));
//}








