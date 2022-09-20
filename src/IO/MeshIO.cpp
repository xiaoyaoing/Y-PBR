#include "MeshIO.hpp"
#include "FileUtils.hpp"

#include "fstream"

namespace  MeshIO{
bool LoadMeshFromFile(const std::string & path,std::vector<Vertex> & vertexs,std::vector<TriangleI>& tris){
        std::ifstream stream(FileUtils::WorkingDir+path);

        if(!stream.is_open())
            return false;

        uint64 numVerts, numTris;
        FileUtils::streamRead(stream, numVerts);
        vertexs.resize(size_t(numVerts));
        FileUtils::streamRead(stream, vertexs);
        FileUtils::streamRead(stream, numTris);
        tris.resize(size_t(numTris));
        FileUtils::streamRead(stream, tris);

        return true;
    }

}