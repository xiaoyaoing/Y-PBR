#include "CurveIO.hpp"
#include "IO/FileUtils.hpp"
#include <fstream>
#include <array>

//Mainly Copyied From Tunsgten
//https://benedikt-bitterli.me/tungsten.html
namespace CurveIO {

    static const std::array < uint8, 8 > FiberMagic{0x80, 0xBF, 0x80, 0x46, 0x49, 0x42, 0x45, 0x52};
    static const uint16 SupportedMajorVersion = 1;
    static const uint16 CurrentMinorVersion = 0;
    static const uint32 SupportedContentType = 0;
    static const size_t FiberValueSize[] = {1, 1, 2, 2, 4, 4, 8, 8, 4, 8};
    enum FiberValueType
    {
        FIBER_INT8   = 0,
        FIBER_UINT8  = 1,
        FIBER_INT16  = 2,
        FIBER_UINT16 = 3,
        FIBER_INT32  = 4,
        FIBER_UINT32 = 5,
        FIBER_INT64  = 6,
        FIBER_UINT64 = 7,
        FIBER_FLOAT  = 8,
        FIBER_DOUBLE = 9,
        FIBER_TYPE_COUNT
    };

    struct FiberAttribute
    {
        uint64 dataLength;
        uint16 attributeFlags;
        uint8 valueType;
        uint8 valuesPerElement;
        std::string attributeName;
        uint64 elementsPresent;

        FiberAttribute(std::ifstream & in)
                : dataLength(FileUtils::streamRead<uint64>(in)),
                  attributeFlags(FileUtils::streamRead<uint16>(in)),
                  valueType(FileUtils::streamRead<uint8>(in)),
                  valuesPerElement(FileUtils::streamRead<uint8>(in)),
                  attributeName(FileUtils::streamRead<std::string>(in)),
                  elementsPresent(0)
        {
            if (valueType < FIBER_TYPE_COUNT)
                elementsPresent = dataLength/(FiberValueSize[valueType]*valuesPerElement);
        }

        ~FiberAttribute(){

        }

        bool matches(const std::string &name, bool perCurve, FiberValueType type, int numValuesPerElement)
        {
            return attributeName == name
                   && ((attributeFlags & 1) != 0) == perCurve
                   && FiberValueType(valueType) == type
                   && valuesPerElement == numValuesPerElement
                   && elementsPresent > 0;
        }

        template<typename T>
        std::unique_ptr<T[]> load(std::ifstream & in, uint64 elementsRequired)
        {
            std::unique_ptr<T[]> result(new T[size_t(elementsRequired)]);
            FileUtils::streamRead(in, result.get(), size_t(elementsPresent));
            // Copy-extend
            for (size_t i = size_t(elementsPresent); i < size_t(elementsRequired); ++i)
                result[i] = result[size_t(elementsPresent - 1)];
            return std::move(result);
        }
    };

    static bool loadFiber(const std::string & path, CurveData & data) {
        std::ifstream in(FileUtils::WorkingDir + path);
        if ( ! in )
            return false;

        if ( FileUtils::streamRead < std::array < uint8, 8>>(in) != FiberMagic )
            return false;

        uint16 versionMajor = FileUtils::streamRead < uint16 >(in);
        uint16 versionMinor = FileUtils::streamRead < uint16 >(in);
        if ( versionMajor != SupportedMajorVersion )
            return false;

        uint32 contentType = FileUtils::streamRead < uint32 >(in);
        if ( contentType != SupportedContentType )
            return false;

        uint64 headerLength = FileUtils::streamRead < uint64 >(in);
        uint64 numVertices = FileUtils::streamRead < uint64 >(in);
        uint64 numCurves = FileUtils::streamRead < uint64 >(in);

        uint64 offset = headerLength;
        while ( true ) {
            in.seekg(size_t(offset), std::ios_base::beg);

            uint64 descriptorLength = FileUtils::streamRead < uint64 >(in);
            if ( descriptorLength == 0 )
                break;

            FiberAttribute attribute(in);

            offset += descriptorLength;
            in.seekg(size_t(offset), std::ios_base::beg);

            if ( data.curveEnds && attribute.matches("num_vertices", true, FIBER_UINT16, 1) ) {
                data.curveEnds->resize(size_t(numCurves));
                auto vertexCounts = attribute.load < uint16 >(in, numCurves);
                for ( size_t i = 0 ; i < size_t(numCurves) ; ++ i )
                    ( * data.curveEnds )[i] = uint32(vertexCounts[i]) + ( i > 0 ? ( * data.curveEnds )[i - 1] : 0 );
            } else if ( data.nodeData && attribute.matches("position", false, FIBER_FLOAT, 3) ) {
                data.nodeData->resize(size_t(numVertices));
                auto pos = attribute.load < vec3 >(in, numVertices);
                for ( size_t i = 0 ; i < size_t(numVertices) ; ++ i )
                    ( * data.nodeData )[i] = vec4(pos[i].x, pos[i].y, pos[i].z, ( * data.nodeData )[i].w);
            } else if ( data.nodeData && attribute.matches("width", false, FIBER_FLOAT, 1) ) {
                data.nodeData->resize(size_t(numVertices));
                auto width = attribute.load < float >(in, numVertices);
                for ( size_t i = 0 ; i < size_t(numVertices) ; ++ i )
                    ( * data.nodeData )[i].w = width[i];
            }

            offset += attribute.dataLength;
        }

//        if ( data.curveEnds && data.nodeData && data.nodeNormal )
//            initializeRandomNormals(data);

        return true;
    }


    bool CurveIO::load(const std::string & path, CurveData & data) {
        auto ext = FileUtils::getFileSuffix(path);
        if(ext == "fiber")
            return loadFiber(path,data);
    }

    bool CurveIO::save(const std::string & path, const CurveIO::CurveData & data) {
        return false;
    }

}
