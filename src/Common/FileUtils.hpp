#include "ostream"

typedef std::shared_ptr<std::ostream> OutputStreamHandle;


class FileUtils {
    static OutputStreamHandle openFileOutputStream(const std::string &p);
};


