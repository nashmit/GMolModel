#include "Simbody.h"
#include "Molmodel.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/QR>

#include <boost/algorithm/string.hpp>

#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cctype>

#include <boost/math/distributions/normal.hpp>

#include "bgeneral.hpp"
#include "SetupReader.hpp"


#ifdef DEBUG_ROBO
#define TRACE(STR) printf("%s", STR);
#else
#define TRACE(STR)
#endif
