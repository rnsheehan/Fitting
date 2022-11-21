#ifndef ATTACH_H
#define ATTACH_H

// Be careful when trying to minimise the librairies that you include
// e.g. iostream includes cmath indirectly at least on MSVS
// https://stackoverflow.com/questions/29338108/is-cmath-or-math-h-really-needed-compiles-without-it
// If you are going to use a library function (or macro/template/whatever), it's up to you to include the correct header.
// Otherwise your program compiling correctly is simply an accident.

#include <cstdlib>
#include <iostream> // cout, cin, cerr
#include <iomanip> // setw, setprecision, time

#include <string>
#include <sstream>
#include <fstream>

// need these for directory manipulation
#include <direct.h>
#include <errno.h>

#include <cmath>
#include <vector>
#include <list>
#include <complex>
#include <iterator>

#include <algorithm>

// Constants
static const double p = (atan(1.0)); // pi / 4
static const double Two_PI = (8.0 * p); // 2 pi
static const double PI = (4.0 * p); // pi
static const double PI_2 = (2.0 * p); // pi / 2
static const double PI_3 = ((4.0 / 3.0) * p); // pi / 3
static const double PI_4 = (p); // pi / 4
static const double PI_5 = ((4.0 / 5.0) * p); // pi / 5
static const double PI_6 = ((2.0 / 3.0) * p); // pi / 6 
static const double EPS = 1.0e-12;
static const double TINY = 1.0e-20;
static const double SQRT_PI = sqrt(PI);

static const int MAX_PATH_LENGTH = 250; // max. length for a directory in Windows OS

static const std::string empty_str = "";
static const std::string dottxt = ".txt";

static const std::complex<double> zero(0.0, 0.0);
static const std::complex<double> eye(0.0, 1.0);
static const std::complex<double> one(1.0, 0.0);

#include "Faddeeva.hh"

#endif
