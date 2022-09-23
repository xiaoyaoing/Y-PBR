
//2022/7/13
#include "Ray.hpp"





Ray::Ray() {

}

Ray::Ray(vec3 start, vec3 direction, Float nt, Float ft)
   : o(start),d(direction),nearT(nt),farT(ft)
{

}
