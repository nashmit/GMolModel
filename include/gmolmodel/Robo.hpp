#include "Simbody.h"
#include "Molmodel.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <boost/algorithm/string.hpp>

#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cctype>

#include <boost/math/distributions/normal.hpp>

#ifndef DEBUG_ROBO
#define DEBUG_ROBO
#endif

#ifdef DEBUG_ROBO
#define TRACE(STR) printf("%s", STR);
#else
#define TRACE(STR)
#endif
