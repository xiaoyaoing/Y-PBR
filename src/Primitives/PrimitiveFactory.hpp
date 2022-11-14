#include "scene.hpp"
#include "Common/Json.hpp"
#include "Common/Transform.hpp"


class Primitive;

namespace PrimitiveFactory{
    void LoadPrimitiveFromJson(const Json & json, Scene & scene);
    void LoadPrimitivesFromJson(const Json & json, Scene & scene);

}