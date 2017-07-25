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

std::vector<int> code(int input);

double decode(std::vector<int> encoded);

std::pair<int, double> R(std::vector<int> &x, double epsilon);

dataSet GenProj(double m, int d);

dataSet PROT_FO(std::vector<double> v, double epsilon, double beta);

std::pair<double, double> PROT_PP_S_Hist_pp(dataSet &data, double epsilon);

int calculateMean(dataSet &data);

int undoNorm(std::pair<std::vector<double>, std::vector<double>> max_min_pair, std::vector<double> &preNorm);
#endif

// PROT_S_HIST
//          PROT_PP_S_Hist_pp
//      PROT_FO
//      AFO
