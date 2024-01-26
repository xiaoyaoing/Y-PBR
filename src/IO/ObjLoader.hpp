#include <fstream>
#include <vector>
#include <unordered_map>
#include "Common/math.hpp"
#include "Primitives/TriangleHelper.hpp"
class ObjLoader {
public:
    ObjLoader(std::ifstream& stream);
    void   loadFile(std::ifstream& stream);
    void   loadLine(const char* line);
    void   loadFace(const char* line);
    uint32 fetchVertex(int32 pos, int32 normal, int32 uv);

public:
    std::vector<vec3>                 _pos, _normal;
    std::vector<vec2>                 _uv;
    std::vector<TriangleI>            _tris;
    std::vector<Vertex>               _verts;
    std::unordered_map<ivec3, uint32> _indices;
    uint32                            _currentMaterial;
    void                              skipWhitespace(const char*& s);

    bool hasPrefix(const char* s, const char* pre);
};