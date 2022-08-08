
//2022/7/13
#include "Ray.hpp"


Ray::Ray(vec3 start, vec3 direction): o(start) , d(direction),
                                      nearT(0), farT(std::numeric_limits<Float>::max()){

}

Ray::Ray() {

}

Ray::Ray(vec3 start, vec3 direction, Float ft, Float nt)
   : o(start),d(direction),nearT(nt),farT(ft)
{

}
