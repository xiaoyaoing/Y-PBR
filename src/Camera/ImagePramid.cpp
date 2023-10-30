#include "ImagePrid.h"

void ImagePrid::saveOutPut(const std::string &fileName, Float scale) {
        auto prefix = FileUtils::getFilePrefix(fileName);
        for (int length = 1; length <= maxPathLength; ++length) {
            for (int t = 1; t <= length + 1; ++t) {

                int s = length + 1 - t;
                if(s!=1 && s!=0)
                    continue;
                int idx = pyramidIndex(s, t);
                char buffer[100];
                sprintf(buffer, "%s-l=%d-c=%d.exr", prefix.c_str(),s,t);
                images[idx].save(std::string(buffer),scale,true);
            }
        }

}