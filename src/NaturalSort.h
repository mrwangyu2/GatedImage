#ifndef natsort_hpp
#define natsort_hpp

#include <string>

/**
* standard C natural string compare
* @param s1 left string
* @param s2 right string
* @return -1 when s1 < s2, 0 when s1 == s2, 1 when s1 > s2
*/
int natstrcmp(const char* s1, const char* s2);

/**
* STL natural less-than string compare
* @param s1 left string
* @param s2 right string
* @return true when natural s1 < s2
*/
bool natstrlt(const char* s1, const char* s2);

/**
* @param s1 left string
* @param s2 right string
* std::string variant of natstrlt.
* @return true when natural s1 < s2
*/
inline bool stlnatstrlt(const std::string& s1, const std::string &s2) {
	return natstrlt(s1.c_str(), s2.c_str());
}

#endif