#include "fstream"
#include "vector"

#pragma once

class  FileUtils{

public:
    static std::string WorkingDir;

    template<typename T>
    static inline void streamRead(std::istream &in, T &dst)
    {
        in.read(reinterpret_cast<char *>(&dst), sizeof(T));
    }

    template<typename T>
    static inline T streamRead(std::istream  &in)
    {
        T t;
        streamRead(in, t);
        return t;
    }

    template<typename T>
    static inline void streamRead(std::istream &in, std::vector<T> &dst)
    {
        in.read(reinterpret_cast<char *>(&dst[0]), dst.size()*sizeof(T));
    }

    template<typename T>
    static inline void streamRead(std::istream &in, T * dst,size_t n)
    {
        in.read(reinterpret_cast<char *>(dst), n*sizeof(T));
    }

    static std::string getFileFullPath(const std::string & path);
    static std::string getFilePath(const std::string & path,const std::string & suffix,bool overwrite = false);
    static std::string getFilePath(const std::string & path,bool overwrite = false);

    static std::string getFileSuffix(const std::string & path);
    static std::string getFilePrefix(const std::string & path);
    static size_t getFileSize(const std::string & path);
};

template<>
inline std::string FileUtils::streamRead<std::string>(std::istream  &in)
{
    std::string s;
    std::getline(in, s, '\0');
    return std::move(s);
}