#include "helpers.h"

#include <iostream>

pair<State, State> extrapolate_state(State s, State d, double t) {
  State s_final(3);
  State d_final(3);

  s_final <<  s[0] + (s[1] * t) + s[2] * pow(t, 2) / 2.0,
              s[1] + s[2] * t,
              s[2];

  d_final <<  d[0] + (d[1] * t) + d[2] * pow(t, 2) / 2.0,
              d[1] + d[2] * t,
              d[2];

  return pair<State, State>(s_final, d_final);
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}


// Evaluate a polynomial.
double polyeval(Coefficients coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

Coefficients JMT(vector< double> start, vector <double> end, double T)
{
    /*
    Calculate the Jerk Minimizing Trajectory that connects the initial state
    to the final state in time T.

    INPUTS

    start - the vehicles start location given as a length three array
        corresponding to initial values of [s, s_dot, s_double_dot]

    end   - the desired end state for vehicle. Like "start" this is a
        length three array.

    T     - The duration, in seconds, over which this maneuver should occur.

    OUTPUT
    an array of length 6, each value corresponding to a coefficent in the polynomial
    s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5
    */

    cout << start[0] << " " << start[1] << " " << start[2] << endl;
    cout << end[0] << " " << end[1] << " " << end[2] << endl;

    MatrixXd tMat(3, 3);
    tMat << pow(T, 3), pow(T, 4), pow(T, 5),
            3*pow(T, 2), 4*pow(T, 3), 5*pow(T, 4),
            6*T, 12*pow(T, 2), 20*pow(T, 3);

    VectorXd sVec(3);
    sVec << end[0] - (start[0] + start[1]*T + 0.5*start[2]*pow(T, 2)),
            end[1] - (start[1] + start[2]*T),
            end[2] - start[2];

    VectorXd alpha = tMat.colPivHouseholderQr().solve(sVec);
    VectorXd coeffs(6);
    coeffs << start[0],start[1],start[2]/2,alpha[0],alpha[1],alpha[2];
    return coeffs;
}
