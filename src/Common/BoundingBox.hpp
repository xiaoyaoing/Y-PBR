#pragma  once
#include "../Common/math.hpp"
#include "../Ray/Ray.hpp"
//#include "glog/logging.h"
class Bounds3 {
public:
    // Bounds3 Public Methods
    Bounds3( ) {
        Float minNum = std::numeric_limits < Float >::lowest();
        Float maxNum = std::numeric_limits < Float >::max();
        pMin = vec3(maxNum);
        pMax = vec3(minNum);
    }

    explicit Bounds3(const vec3 & p) : pMin(p), pMax(p) {}

    Bounds3(const vec3 & p1, const vec3 & p2)
            : pMin(min(p1,p2)),
              pMax(max(p1,p2)) {}

    const vec3 & operator [](int i) const {
        assert(i==0 || i==1);
        return i==0?pMin:pMax;
    }

    vec3 & operator [](int i){
        assert(i==0 || i==1);
        return i==0?pMin:pMax;
    }

    bool operator ==(const Bounds3  & b) const {
        return b.pMin == pMin && b.pMax == pMax;
    }

    bool operator !=(const Bounds3  & b) const {
        return b.pMin != pMin || b.pMax != pMax;
    }

    vec3 Corner(int corner) const {
//        DCHECK(corner >= 0 && corner < 8);
        return {( * this )[( corner & 1 )].x,
                            ( * this )[( corner & 2 ) ? 1 : 0].y,
                            ( * this )[( corner & 4 ) ? 1 : 0].z};
    }

    vec3 Diagonal( ) const { return pMax - pMin; }

    Float SurfaceArea( ) const {
        vec3 d = Diagonal();
        return 2 * ( d.x * d.y + d.x * d.z + d.y * d.z );
    }

    Float Volume( ) const {
        vec3 d= Diagonal();
        return d.x * d.y * d.z;
    }

    int MaximumExtent( ) const {
        vec3 d = Diagonal();
        if ( d.x > d.y && d.x > d.z )
            return 0;
        else if ( d.y > d.z )
            return 1;
        else
            return 2;
    }

//    vec3 Lerp(const vec3 & t) const {
//        return lerp(pMin,pMax,t);
//    }

    vec3 Offset(const vec3  & p) const {
        vec3 o = p - pMin;
        if ( pMax.x > pMin.x ) o.x /= pMax.x - pMin.x;
        if ( pMax.y > pMin.y ) o.y /= pMax.y - pMin.y;
        if ( pMax.z > pMin.z ) o.z /= pMax.z - pMin.z;
        return o;
    }

//    void BoundingSphere(vec3 * center, Float * radius) const {
//        * center = ( pMin + pMax ) / 2;
//        * radius = Inside(* center, * this) ? Distance(* center, pMax) : 0;
//    }

//    template < typename U >
//    explicit operator Bounds3 < U >( ) const {
//        return Bounds3 < U >((Point3 <U>) pMin, (Point3 <U>) pMax);
//    }
    bool Contains(const vec3 & p) const {
        return (p.x >= pMin.x && p.x <= pMax.x && p.y >= pMin.y &&
                p.y <= pMax.y && p.z >= pMin.z && p.z <= pMax.z);
    }

    void BoundingSphere(vec3 *center, Float *radius) const {
        *center = (pMin + pMax) /2.f ;
        *radius = Contains(*center) ?distance(*center, pMax) : 0;
    }

    Float BoundingRadius(){
        return distance(pMin, pMax) * 0.25;
    }

     bool IntersectP(const Ray & ray, const vec3 & invDir,
                           const int dirIsNeg[3]) const;

//    friend std::ostream & operator <<(std::ostream & os, const Bounds3 < T > & b) {
//        os << "[ " << b.pMin << " - " << b.pMax << " ]";
//        return os;
//    };

    vec3 pMin,pMax;
};

inline  Bounds3  Union(const Bounds3 & a,const Bounds3 & b){
    return Bounds3(min(a.pMin,b.pMin),max(a.pMax,b.pMax));
}

inline  Bounds3  Union(const Bounds3 & a,const vec3 & v){
    return Bounds3(min(a.pMin,v),max(a.pMax,v));
}

inline  Bounds3  Expand(const Bounds3 & a,Float length){
    return Bounds3(a.pMin-length,a.pMax + length);
}

inline  bool Inside(const Bounds3 & b,const vec3 & p){
    return (p.x >= b.pMin.x && p.x <= b.pMax.x && p.y >= b.pMin.y &&
            p.y <= b.pMax.y && p.z >= b.pMin.z && p.z <= b.pMax.z);
}


