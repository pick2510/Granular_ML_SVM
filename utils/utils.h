#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

template <class Container>
void split(const std::string& str, Container& cont)
{
    std::istringstream iss(str);
    std::copy(std::istream_iterator<std::string>(iss),
         std::istream_iterator<std::string>(),
         std::back_inserter(cont));
}

#endif