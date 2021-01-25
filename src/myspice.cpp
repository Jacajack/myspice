#include "legacy.hpp"
#include "extended.hpp"

/**
	\file myspice.cpp
	\brief Główny plik programu
	\author Jacek Wieczorek
*/

int main(int argc, char *argv[])
{
	#ifdef EXTENDED_MODE
		return main_extended(argc, argv);
	#else
		return main_legacy(argc, argv);
	#endif
}