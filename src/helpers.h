#ifndef HELPERS_H
#define HELPERS_H
#pragma once

#include <math.h>
#include <vector>
#include <tuple>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using trajectory = vector<double>; // length 6
using CostFunctionPointer = double(*)(trajectory); // type alias
using CostPair = pair<CostFunctionPointer, double>;

// Coefficients_S, Coefficients_D, duration_time
using Coefficients = Eigen::VectorXd; // length 6;
using Path = tuple<Coefficients, Coefficients, double>;

// x, x_dot, x_dot_dot
using State = Eigen::VectorXd;

// Evaluate a polynomial.
double polyeval(Coefficients coeffs, double x);

pair<State, State> extrapolate_state(State s, State d, double t);

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order);

Coefficients JMT(vector<double> start, vector <double> end, double T);
Path PTG(State start_s, State start_d);

//   perturb_goal
//   calculate_cost
// }

#endif /* HELPERS_H */
