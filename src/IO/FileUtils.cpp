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


