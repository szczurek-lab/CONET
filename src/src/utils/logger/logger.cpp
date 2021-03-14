#include "logger.h"

std::ofstream LOG_FILE{ "cnv_logs", std::ios_base::app };

void print()
{
	LOG_FILE << "\n";
}

void printSTDOUT()
{
	std::cout << "\n";
}

void printErr()
{
	LOG_FILE << "\n";
}
