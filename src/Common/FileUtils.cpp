#include "FileUtils.hpp"

//
//OutputStreamHandle FileUtils::openFileOutputStream(const std::string & p) {
//
//    std::shared_ptr<std::ostream> out(new std::ofstream(p,
//                                                        std::ios_base::out | std::ios_base::binary),
//                                      [](std::ostream *stream){ finalizeStream(stream); });
//
//    OutputStreamHandle out = openFileOutputStream(p);
//
//
//
//    return out;
//
//}