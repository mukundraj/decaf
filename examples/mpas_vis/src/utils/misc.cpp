#include <iostream>
#include <string>
#include <sstream>
#include "misc.h"

std::string itos(int i) // convert int to string
{
	std::stringstream s;
	s << i;
	return s.str();

}


ghost_req::~ghost_req(){}

std::string split_filename(std::string str){
	std::size_t found = str.find_last_of("/\\");
	// std::string fullpath = str;
	return str.substr(found+1);
}