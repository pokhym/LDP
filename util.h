#ifndef __UTIL_H
#define __UTIL_H

#include "data.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>    // exp
#include <random>   // bernoulli



int parseData(char const* filename, dataSet &data);

int normalizeNeg1toPos1(dataSet &data, std::vector<double> outlier);

std::vector<double> tuplePerturbation(dataSet &data, double epsilon);

int calculateMean(dataSet &data);

#endif
