#include "fstream"
#include "vector"



class  FileUtils{

public:
    static std::string WorkingDir;

    template<typename T>
    static inline void streamRead(std::istream &in, T &dst)
    {
        in.read(reinterpret_cast<char *>(&dst), sizeof(T));
    }

    template<typename T>
    static inline void streamRead(std::istream &in, std::vector<T> &dst)
    {
        in.read(reinterpret_cast<char *>(&dst[0]), dst.size()*sizeof(T));
    }

    static std::string getFilePath(const std::string & path,const std::string & suffix,bool overwrite);
};