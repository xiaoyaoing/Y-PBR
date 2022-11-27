#pragma  once

#include "SampleRecords/SurfaceScatterEvent.hpp"
#include "Colors/Spectrum.hpp"
#include "BsdfTypes.hpp"
#include "Common/Texture.hpp"

#include "Common/Json.hpp"
#include <spdlog/spdlog.h>


// BSDF Inline Functions
inline Float CosTheta(const vec3 & w) { return w.z; }

inline Float Cos2Theta(const vec3 & w) { return w.z * w.z; }

inline Float AbsCosTheta(const vec3 & w) { return std::abs(w.z); }

inline Float Sin2Theta(const vec3 & w) {
    return std::max((Float) 0, (Float) 1 - Cos2Theta(w));
}

inline Float SinTheta(const vec3 & w) { return std::sqrt(Sin2Theta(w)); }

inline Float TanTheta(const vec3 & w) { return SinTheta(w) / CosTheta(w); }

inline Float Tan2Theta(const vec3 & w) {
    return Sin2Theta(w) / Cos2Theta(w);
}

inline Float CosPhi(const vec3 & w) {
    Float sinTheta = SinTheta(w);
    return ( sinTheta == 0 ) ? 1 : std::clamp(w.x / sinTheta, - 1.f, 1.f);
}

inline Float SinPhi(const vec3 & w) {
    Float sinTheta = SinTheta(w);
    return ( sinTheta == 0 ) ? 0 : std::clamp(w.y / sinTheta, - 1.f, 1.f);
}

inline Float Cos2Phi(const vec3 & w) { return CosPhi(w) * CosPhi(w); }

inline Float Sin2Phi(const vec3 & w) { return SinPhi(w) * SinPhi(w); }


inline vec3 Reflect(const vec3 & wo, const vec3 & n) {
    return - wo + 2 * dot(wo, n) * n;
}

inline bool SameHemisphere(const vec3 & w, const vec3 & wp) {
    return w.z * wp.z > 0;
}


//class BSDF;
////
//class Bsdf{
//
//public:
//    Spectrum f(const vec3 &wo , const vec3 &wi,
//               BXDFType flags = BSDF_ALL) const ;
//
//    ///sample a new dir
//    ///wi wo in local space
//    virtual Spectrum sampleF(const vec3 &wo, vec3 *wi, const vec2 &u,
//                             Float *pdf, BXDFType type,BXDFType *sampledType) const ;
//
//    int NumComponents(BXDFType flags = BSDF_ALL) const;
//
//    void Add(BXDF *b) {
//        BXDFs[nBXDFs++] = b;
//    }
//    void LogInfo();
//    std::string name;
//
////    Float eta(){
////        return BXDF[]
////    }
//private:
//    BXDF *BXDFs[8];
//    int nBXDFs = 0;
//    //for debug
//};

class BSDF {
public:

    BSDF(BXDFType type) : m_type(type), m_albedo(nullptr), m_bumpMap(nullptr) {};




    virtual Float Pdf(const SurfaceScatterEvent & event) const = 0;


    inline Spectrum sampleF(SurfaceScatterEvent & event, const vec2 & u, bool adjoint) const  {

        Spectrum fResult = sampleF(event, u);

        if(adjoint ) return fResult;

        fResult *= sqr(eta(event));
        return fResult;
    }

    inline Spectrum f(const SurfaceScatterEvent & event, bool adjoint) const {
        Spectrum fResult = f(event);
        fResult *= sqr(eta(event));
        return fResult;
    }


    virtual void LogInfo( ) const = 0;

    bool MatchesFlags(BXDFType t) const {
        return ( m_type & t ) == m_type;
    }

    bool hasFlag(BXDFType t) const {
        return ( m_type & t ) == t;
    }

    virtual Float eta(const SurfaceScatterEvent & event) const { return 1; }

    void setAlbedo(const std::shared_ptr < Texture < Spectrum>> albedo) {
        m_albedo = std::move(albedo);
    }


    void setBumpMap(const std::shared_ptr < Texture < Float>> & bumpMap) { m_bumpMap = bumpMap; }

    std::string name;
    std::shared_ptr < Texture < Spectrum>> m_albedo = nullptr;
    int temp;
protected:
    virtual Spectrum sampleF(SurfaceScatterEvent & event, const vec2 & u) const = 0;
    virtual Spectrum f(const SurfaceScatterEvent & event) const = 0;


    std::shared_ptr < Texture < Float>> m_bumpMap = nullptr;
    BXDFType m_type;
};


class Mirror : public BSDF {

};

class LambertainR : public BSDF {

public:
    LambertainR( ) : BSDF(BXDFType(BSDF_DIFFUSE | BSDF_REFLECTION)) {}

    virtual Spectrum f(const SurfaceScatterEvent & event) const override;

    Float Pdf(const SurfaceScatterEvent & event) const override;

    virtual void LogInfo( ) const;


    virtual Spectrum sampleF(SurfaceScatterEvent & event, const vec2 & u) const override;

private:
    bool useCosineSample;
};

class LambertainT : public BSDF {
public:

    virtual Spectrum f(const SurfaceScatterEvent & event) const override;

    virtual void LogInfo( ) const override;

    LambertainT( ) : BSDF(BXDFType(BSDF_DIFFUSE | BSDF_TRANSMISSION)) {}


    virtual Spectrum sampleF(SurfaceScatterEvent & event, const vec2 & u) const override;

private:
};

class SpecularR : public BSDF {
public:


    virtual Spectrum f(const SurfaceScatterEvent & event) const override;

    virtual void LogInfo( ) const override;

    SpecularR( ) : BSDF(BXDFType(BSDF_REFLECTION | BSDF_SPECULAR)) {}

    virtual Spectrum sampleF(SurfaceScatterEvent & event, const vec2 & u) const override;

    // avoid calculate specular pdf directly
    Float Pdf(const SurfaceScatterEvent & event) const override;
};


class OrenNayar : public BSDF {
public:
    OrenNayar(const Spectrum & R, Float sigma);

private:
    const Spectrum R;
    Float A, B;
};




