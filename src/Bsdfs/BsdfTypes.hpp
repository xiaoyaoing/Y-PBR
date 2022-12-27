// BSDF Declarations
#pragma  once

enum BXDFType {
    BSDF_REFLECTION = 1 << 0,
    BSDF_TRANSMISSION = 1 << 1,
    BSDF_DIFFUSE = 1 << 2,
    BSDF_GLOSSY = 1 << 3,
    BSDF_SPECULAR = 1 << 4,
    BSDF_FORWARD = 1<<5,
    BSDF_ALL = BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR | BSDF_REFLECTION |
               BSDF_TRANSMISSION,
   // BSDF_NO_SPECULAR = ~ (BSDF_SPECULAR | BSDF_FORWARD)
    BSDF_NO_SPECULAR = BSDF_DIFFUSE | BSDF_GLOSSY
};

inline bool isSpecualr(BXDFType type){
    return BXDFType(type &  BSDF_SPECULAR) ==  BSDF_SPECULAR;
}