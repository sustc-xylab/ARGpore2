// Copyright 2017 Martin C. Frith

#ifndef GETOPT_UTIL_HH
#define GETOPT_UTIL_HH

#include <unistd.h>

inline void resetGetopt() {  // XXX fragile voodoo
#ifdef __GLIBC__
  optind = 0;
#else
  optind = 1;
  //optreset = 1;  // XXX ???
#endif
}

#endif
