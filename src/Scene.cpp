#include "scene.hpp"
#include "IO/MeshIO.hpp"
#include "IO/FileUtils.hpp"

#include "Bsdfs/Reflection.hpp"
#include "Bsdfs/BsdfFactory.hpp"

#include "Primitives/PrimitiveFactory.hpp"
#include "Primitives/TriangleMesh.hpp"
#include "Primitives/Sphere.hpp"
#include "Primitives/Quad.hpp"
#include "Primitives/Cube.hpp"

#include "Lights/Infinte.hpp"
#include "Lights/Distant.hpp"
#include "Lights/SkyDome.hpp"
#include "Lights/InfiniteSphereCap.h"

#include "Texture/TextureFactory.hpp"
#include "Common/Transform.hpp"

#include <spdlog/spdlog.h>
#include <iostream>

static int occludedCount = 0;
static int occludedCount1 = 0;


Scene::Scene(const Json sceneJson) : options(RenderOptions(sceneJson.at("renderer"))) {
    bsdfs = BSDFFactory::LoadBsdfsFromJson(sceneJson.at("bsdfs"));
    PrimitiveFactory::LoadPrimitivesFromJson(sceneJson.at("primitives"), * this);

    spdlog::info("{} Primitives", primitives.size());
    spdlog::info("{} lights", lights.size());
    spdlog::info("{} Bsdfs", bsdfs.size());

    _useBVH = getOptional(sceneJson["renderer"], "scene_bvh", false);

}


void Scene::handleAddLight(const Json & p, int l, int r) {
    if ( contains(p, "emission") ) {
        std::shared_ptr < Texture < Spectrum>> emission = TextureFactory::LoadTexture < Spectrum >(p, "emission",
                                                                                                   Spectrum(1));
        for ( size_t i = l ; i < r ; i ++ ) {
            auto light = std::make_shared < AreaLight >(this->primitives[i], emission);
            this->primitives[i]->areaLight = light;
            lights.push_back(light);
        }
    } else if ( contains(p, "power") ) {
        std::shared_ptr < Texture < Spectrum>> power = TextureFactory::LoadTexture < Spectrum >(p, "power",
                                                                                                Spectrum(1));
        Float totalScale = 0;
        for ( size_t i = l ; i < r ; i ++ ) {
            totalScale += 1 / primitives[i]->powerToRadianceScale();
        }
        totalScale = 1 / totalScale;
        power->setScale(totalScale);

        for ( size_t i = l ; i < r ; i ++ ) {
            auto light = std::make_shared < AreaLight >(this->primitives[i], power);
            this->primitives[i]->areaLight = light;
            lights.push_back(light);
        }
    }

    return;
}


std::optional < Intersection > Scene::intersect(const Ray & ray) const {
    if ( _useBVH ) {
        RTCRayHit rayHit;
        EmbreeUtils::convertRay(& ray, & rayHit);
        //return bvh->intersect(ray);
        RTCIntersectContext context;
        rtcInitIntersectContext(& context);
        rtcIntersect1(_scene, & context, & rayHit);
        if ( rayHit.hit.geomID == RTC_INVALID_GEOMETRY_ID ) {
            return std::nullopt;
        }
        Intersection * its = EmbreeUtils::RTCRayHit1_(& rayHit)->its;
        its->w = ray.d;
        return {* its};
    }

    std::optional < Intersection > minIntersection;
    Ray _ray(ray);
    for ( auto primitive: primitives ) {

        auto its = primitive->intersect(_ray);
        if ( its.has_value() ) {
            minIntersection = its;
        }
    }
    if ( minIntersection.has_value() ) {
        minIntersection->w = ray.d;
    }
    return minIntersection;
}

bool Scene::intersectP(const Ray & ray) const {
    if ( _useBVH ) {
        RTCRay shadowRay;
        EmbreeUtils::convertRay(& ray, & shadowRay);
        RTCIntersectContext context;
        rtcInitIntersectContext(& context);
        rtcOccluded1(_scene, & context, & shadowRay);
        return shadowRay.tfar != ray.farT;
    } else {
        Ray _ray(ray);
        for ( auto primitive: primitives ) {
            if ( primitive->occluded(ray) )
                return true;
        }
    }

    return false;
}

void Scene::logDebugInfo( ) {
    spdlog::info("Occluded {0} {1}", occludedCount, occludedCount1);
//    static_cast<InfinteSphere *>(lights[0].get())->logDebugInfo();
}


void Scene::build( ) {


    if ( _useBVH ) {
        _scene = rtcNewScene(EmbreeUtils::getDevice());
        for ( const auto & primitve: primitives ) {
            RTCGeometry geometry = primitve->initRTC();
            auto geom = rtcNewGeometry(EmbreeUtils::getDevice(), RTC_GEOMETRY_TYPE_USER);
            rtcEnableGeometry(geom);
            rtcSetGeometryUserPrimitiveCount(geom, 1);
            rtcSetGeometryUserData(geom, primitve.get());
            rtcSetGeometryBoundsFunction(geom, & EmbreeUtils::instanceBoundsFunc, nullptr);
            rtcSetGeometryIntersectFunction(geom, & EmbreeUtils::instanceIntersectFunc);
            rtcSetGeometryOccludedFunction(geom, & EmbreeUtils::instanceOccludedFunc);
            rtcCommitGeometry(geom);
            int res = rtcAttachGeometry(_scene, geom);
            spdlog::info(res);
        }
        rtcCommitScene(_scene);
        //bvh = std::make_unique<BVHAccel>(primitives);  // I will use embree now ,the bvh I implmented is too slow!
    }

    for ( auto p: primitives ) {
        worldBound = Union(worldBound, p->BB());
    }
    for ( auto light: lights )
        light->Preprocess(* this);

}
