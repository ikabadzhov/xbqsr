#ifndef GENOMICPIPELINE_QNAMESTRINGMAPPER_H
#define GENOMICPIPELINE_QNAMESTRINGMAPPER_H

#include <vector>
#include <string>
#include <algorithm>

class QnameStringMapper {
private:
    size_t REALLOCATION = 5000000;
    ssize_t nextPosition{};
    std::vector<std::string> sortedStrings;
    std::vector<ssize_t> firstStringsIdxs;

public:
    /*static QnameStringMapper& getInstance() {
        static QnameStringMapper instance;
        return instance;
    }*/
    QnameStringMapper() = default;

    ssize_t addString(const std::string& str);
    void resize(size_t sz);
    void finalize(std::vector<int32_t> &qname);
    const std::string& GetAt(int32_t id);

    static const ssize_t EMPTY_STRING = -1;
};



#endif //GENOMICPIPELINE_QNAMESTRINGMAPPER_H
