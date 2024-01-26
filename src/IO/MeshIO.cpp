#include "MeshIO.hpp"
#include "FileUtils.hpp"
#include "ObjLoader.hpp"

#include <fstream>

///file format wo3 from Tungsten
namespace MeshIO {
    bool LoadWo3(std::ifstream& stream, std::vector<Vertex>& vertexs, std::vector<TriangleI>& tris) {
        uint64 numVerts, numTris;
        FileUtils::streamRead(stream, numVerts);
        vertexs.resize(size_t(numVerts));
        FileUtils::streamRead(stream, vertexs);
        FileUtils::streamRead(stream, numTris);
        tris.resize(size_t(numTris));
        FileUtils::streamRead(stream, tris);

        return true;
    }

    bool LoadObj(std::ifstream& stream, std::vector<Vertex>& vertexs, std::vector<TriangleI>& tris) {
        ObjLoader objLoader(stream);
        vertexs = std::move(objLoader._verts);
        tris    = std::move(objLoader._tris);
        return true;
    }
    bool LoadMeshFromFile(const std::string& path, std::vector<Vertex>& vertexs, std::vector<TriangleI>& tris) {
        std::ifstream stream(FileUtils::WorkingDir + path, std::ios::binary);
        if (!stream.is_open())
            return false;
        if (FileUtils::getFileSuffix(path) == "wo3")
            return LoadWo3(stream, vertexs, tris);
        if (FileUtils::getFileSuffix(path) == "obj")
            return LoadObj(stream, vertexs, tris);
        return false;
    }

}