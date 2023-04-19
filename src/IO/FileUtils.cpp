#include <sys/stat.h>
#include "FileUtils.hpp"

std::string FileUtils::WorkingDir = "";

inline bool fileExists (const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

template<typename ... Args>
std::string stringFormat( const char * format, Args ... args )
{
    int size_s = std::snprintf( nullptr, 0, format, args ... ) + 1;
    if( size_s <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
    auto size = static_cast<size_t>( size_s );
    std::unique_ptr<char[]> buf( new char[ size ] );
    std::snprintf( buf.get(), size, format, args ... );
    return std::string( buf.get(), buf.get() + size - 1 );
}


std::string FileUtils::getFilePath(const std::string & path,const std::string & suffix,bool overwrite){
    std::string destPath = path + "." +suffix;
    if(overwrite)
        return std::move(destPath);
    int count = 1;
    while( fileExists(destPath)){
        destPath = stringFormat("%s%d.%s",path.c_str(),count++,suffix.c_str());
    }
    return std::move(destPath);
}

size_t FileUtils::getFileSize(const std::string & path) {
    struct stat statbuf;
    stat(path.c_str(), &statbuf);
    size_t filesize = statbuf.st_size;
    return filesize;
}

std::string FileUtils::getFileSuffix(const std::string & path) {
    int dotIdx = path.find(".");
    if(dotIdx == -1)
        return "";
    return path.substr(dotIdx+1,path.size());
}

std::string FileUtils::getFilePrefix(const std::string & path){
    int dotIdx = path.find(".");
    if(dotIdx == -1)
        return path;
    return path.substr(0,dotIdx);
}

std::string FileUtils::getFilePath(const std::string &path, bool overwrite) {
    return getFilePath(getFilePrefix(path), getFileSuffix(path),overwrite);
}

std::string FileUtils::getFileFullPath(const std::string &path) {
    return WorkingDir+path;
}


