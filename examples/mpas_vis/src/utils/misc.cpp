#include <iostream>
#include <string>
#include <sstream>

std::string itos(int i) // convert int to string
{
	std::stringstream s;
	s << i;
	return s.str();

}
