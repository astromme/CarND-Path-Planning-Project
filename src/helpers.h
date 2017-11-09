#ifndef HELPERS_H
#define HELPERS_H
#pragma once

#include <math.h>
#include <vector>
#include <tuple>
#include "json.hpp"
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

enum DriveState {
  Straight = 0,
  MergeLeft = 1,
  MergeRight = 2,
};

using json = nlohmann::json;

const double dt = 0.02; // 20 ms

using trajectory = vector<double>; // length 6
using CostFunctionPointer = double(*)(trajectory); // type alias
using CostPair = pair<CostFunctionPointer, double>;

// Coefficients_S, Coefficients_D, duration_time
using Coefficients = Eigen::VectorXd; // usually length 6;
using Path = tuple<Coefficients, Coefficients, double>;

// x, x_dot, x_dot_dot
using State = vector<double>;
using DState = State;
using SState = State;
using SDState = pair<SState, DState>;
using Goal = tuple<SState, DState, double>;

// Evaluate a polynomial.
double polyeval(Coefficients coeffs, double x);
Coefficients differentiate(Coefficients coeff);

SDState extrapolate_state(State s, State d, double t);

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order);

Coefficients JMT(vector<double> start, vector <double> end, double T);
Path PTG(DriveState state, State start_s, State start_d, json sensor_fusion, int T);
pair<State, State> perturb_goal(State goal_s, State goal_d);
double calculate_cost(Coefficients s_coeff, Coefficients d_coeff, double T, json sensor_fusion);

#endif /* HELPERS_H */
