#include "PrimitiveFactory.hpp"

#include "Quad.hpp"
#include "Sphere.hpp"
#include "Cube.hpp"
#include "Curve.hpp"
#include "Disk.hpp"

#include "Lights/InfiniteSphere.hpp"
#include "Lights/InfiniteSphereCap.h"
#include "Lights/SkyDome.hpp"
#include "Lights/Distant.hpp"

#include "Texture/TextureFactory.hpp"
#include "TriangleMesh.hpp"

#include <spdlog/spdlog.h>

static std::shared_ptr < Texture < Spectrum>> getEmission(const Json & json, std::shared_ptr < const Primitive > prim) {
    if ( json.contains("emission") )
        return TextureFactory::LoadTexture < Spectrum >(json, "emission");
    if ( json.contains("power") ) {
        auto powerTexture = TextureFactory::LoadTexture < Spectrum >(json, "power");
        powerTexture->setScale(prim->powerToRadianceScale());
        return powerTexture;
    }
    return nullptr;
}

namespace PrimitiveFactory {
    using LoadPrimitve = std::function < void(const Json &, Scene &) >;

    template < class T >
    std::shared_ptr < T > LoadSimpleHelper(const Json & json) { return std::make_shared < T >(json); }

    std::unordered_map < std::string, std::function < std::shared_ptr < Primitive >(const Json &)>> loadSimpleMap = {
            {"quad",   LoadSimpleHelper < Quad >},
            {"sphere", LoadSimpleHelper < Sphere >},
            {"cube",   LoadSimpleHelper < Cube >},
            {"mesh", LoadSimpleHelper <TriangleMesh>}
    };

    std::unordered_map < std::string, std::function < std::shared_ptr < Light >(const Json &)>> loadInfiniteMap = {
            {"infinite_sphere",     LoadSimpleHelper < InfinteSphere >},
            {"infinite_sphere_cap", LoadSimpleHelper < InfinteSphereCap >},
            {"skydome",             LoadSimpleHelper < SkyDome >},
            {"distant",             LoadSimpleHelper < DistantLight >}
    };


    void LoadSimpleFromJson(const Json & json, Scene & scene) {
        auto type = getOptional(json, "type", std::string("sphere"));
        auto prim = loadSimpleMap[type](json);
        prim->load(json, scene);
        prim->transform();
     //   if(type!="mesh")
        scene.AddPrimitive(prim);
        auto emission = getEmission(json, prim);
        if ( emission ) {
            auto light = std::make_shared < AreaLight >(prim, emission);
            prim->areaLight = light;
            scene.AddLight(light);
        }
    }

//    void LoadMeshFromJson(const Json & json, Scene & scene) {
//        const std::string materialName = json["bsdf"];
//
//        std::shared_ptr < TriangleMesh > prim = std::make_shared < TriangleMesh >();
//        prim->build(json, scene);
//        if ( materialName == "Floor" || materialName == "Walls" )
//            scene.AddPrimitive(prim);
//        auto emission = getEmission(json, prim);
//        if ( emission ) {
//            auto light = std::make_shared < AreaLight >(prim, emission);
//            prim->areaLight = light;
//            scene.AddLight(light);
//        }
//    }

    void LoadCurvesFromJson(const Json & json, Scene & scene) {
        std::shared_ptr < Curve > curve = std::make_shared < Curve >(json, scene);
        curve->load(json, scene);
        scene.AddPrimitive(curve);
    }

    void LoadInfiniteFromJson(const Json & json, Scene & scene) {
        auto type = json["type"];
        auto light = loadInfiniteMap[type](json);
        scene.AddLight(light);
    }

    std::unordered_map < std::string, LoadPrimitve > loadPrimMap = {
            {"quad",                LoadSimpleFromJson},
            {"sphere",              LoadSimpleFromJson},
            {"cube",                LoadSimpleFromJson},
            {"mesh",                LoadSimpleFromJson},
            {"infinite_sphere",     LoadInfiniteFromJson},
            {"infinite_sphere_cap", LoadInfiniteFromJson},
            {"skydome",             LoadInfiniteFromJson},
            {"distant",             LoadInfiniteFromJson},
            {"curves",              LoadCurvesFromJson}
    };

    void LoadPrimitiveFromJson(const Json & json, Scene & scene) {
        if ( ! json.contains("type") ) {
            spdlog::error("Primitive without type");
            return;
        }
        std::string type = json["type"];
        if ( ! loadPrimMap.contains(type) ) {
            spdlog::error("Type {0} not supported", type);
            return;
        }
        loadPrimMap[type](json, scene);
    }

    void LoadPrimitivesFromJson(const Json & json, Scene & scene) {
        for ( const Json & subJson: json )
            LoadPrimitiveFromJson(subJson, scene);
    }

}
