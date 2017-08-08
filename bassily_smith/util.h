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

std::vector<double> tuplePerturbation(dataSet &data, double epsilon);

void gen_rand_mtx(boolSet &rand_mtx, int m, int k);

std::vector<double> bassily_smith(dataSet &data, double epsilon, double beta, int k);

unsigned long long nCr(long long n, long long r);

int calculateMean(dataSet &data);

int undoNorm(std::pair<std::vector<double>, std::vector<double>> max_min_pair, std::vector<double> &preNorm);
#endif
