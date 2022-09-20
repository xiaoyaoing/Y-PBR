#include "Common/math.hpp"
#include "Primitives/TriangleHelper.hpp"
#include "vector"
namespace  MeshIO{
    bool LoadMeshFromFile(const std::string &  path,std::vector<Vertex> &  vertexs,std::vector<TriangleI> & tris);
}