#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>

/**
 * Set this to <code>true</code> if debug information should be printed during
 * sampling\n Logging is by default directed to the standard output but its
 * quite easy to modify this code and redirect the output.
 */
constexpr bool DEBUG = false;

void print();

template <class A0, class... Args> void print(A0 a0, Args... args) {
  std::cout << a0;
  print(args...);
}

template <class... Args> void log(Args... args) { return print(args...); }

template <class... Args> void log_err(Args... args) {
  return print("Error: ", args...);
}

template <class... Args> void log_debug(Args... args) {
  if (DEBUG) {
    print("Debug: ", args...);
  }
}

#endif // !LOGGER_H
