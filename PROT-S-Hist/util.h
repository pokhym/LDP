#ifndef __UTIL_H
#define __UTIL_H

#include "data.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>    // exp
#include <random>   // bernoulli



int parseData(char const* filename, dataSet &data, int init_m, int init_n);

std::pair<std::vector<double>, std::vector<double>> normalizeNeg1toPos1(dataSet &data, std::vector<double> outlier);

std::vector<double> e_basic_normalizer(std::vector<double> x, double epsilon);

std::vector<double> tuplePerturbation(dataSet &data, double epsilon);

unsigned long long nCr(long long n, long long r);

int calculateMean(dataSet &data);

int undoNorm(std::pair<std::vector<double>, std::vector<double>> max_min_pair, std::vector<double> &preNorm);
#endif
