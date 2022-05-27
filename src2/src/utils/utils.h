#ifndef UTILS_H
#define UTILS_H

#include <vector> 
#include <algorithm>

namespace Utils {
    template<class T> std::vector<T> flatten(std::vector<std::vector<T>> vec) {
        std::vector<T> flattened; 
        std::for_each(vec.begin(), vec.end(), [&flattened](std::vector<T> &v){
            std::for_each(v.begin(), v.end(), [&flattened](T e) {flattened.push_back(e);});
        });
        return flattened;
    }

}

#endif 