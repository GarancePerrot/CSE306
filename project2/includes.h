#ifndef INCLUDES_H  
#define INCLUDES_H


#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>
#include <random>
#include <ctime>
#include <cstring>
#include <omp.h>
#include <list>
#include <iostream>
#include <sstream>
#include <filesystem>
#include <fstream>
#include <cstdio> 
#include "lbfgs.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

static std::default_random_engine engine(200) ; 
static std::uniform_real_distribution<double> uniform(0.,double(RAND_MAX)-1) ;

#endif 
