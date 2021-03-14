#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <fstream>

extern std::ofstream LOG_FILE;
/**
* Set this to <code>true</code> if debug information should be printed during sampling\n 
* Logging is by default directed to the standard output but its quite easy to modify this code and 
* redirect the output.
*/
constexpr bool DEBUG = false;
/**
* Set this to <code>true</code> if integrity of data structures should be tested during sampling
*/
constexpr bool TEST = false;


void print();

void printErr();

void printSTDOUT();

template <class A0, class ...Args> void print(A0 a0, Args ...args)
{
	LOG_FILE << a0;
	print(args...);
}

template <class A0, class ...Args> void printSTDOUT(A0 a0, Args ...args)
{
	std::cout << a0;
	printSTDOUT(args...);
}

template <class A0, class ...Args> void printErr(A0 a0, Args ...args)
{
	LOG_FILE << a0;
	printErr(args...);
}

template <class ...Args> void log(Args ...args)
{
	return print(args...);
}

template <class ...Args> void logErr(Args ...args)
{
	return printErr("Error: ", args...);
}

template <class ...Args> void logDebug(Args ...args)
{
	return printSTDOUT("Error: ", args...);
}

#endif // !LOGGER_H
