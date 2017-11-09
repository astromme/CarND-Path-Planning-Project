#include "helpers.h"

#include <iostream>
#include <random>


pair<State, State> extrapolate_state(State s, State d, double t) {
  State s_final = { s[0] + (s[1] * t) + s[2] * pow(t, 2) / 2.0,
                    s[1] + s[2] * t,
                    s[2] };

  State d_final = { d[0] + (d[1] * t) + d[2] * pow(t, 2) / 2.0,
                    d[1] + d[2] * t,
                    d[2] };

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

    // cout << start[0] << " " << start[1] << " " << start[2] << endl;
    // cout << end[0] << " " << end[1] << " " << end[2] << endl;

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
    // coeffs << start[0],start[1],start[2]/2,alpha[0],alpha[1],alpha[2];
    return coeffs;
}

double find_lane(double d) {
  return floor(d / 4.0);
}

bool in_same_lane(double d1, double d2) {
  return fabs(d1 - d2) <= 2;
}

pair<State, State> perturb_goal(State goal_s, State goal_d) {
  /*
  Returns a "perturbed" version of the goal.
  */
  static std::default_random_engine gen;
  static std::normal_distribution<double> s0(0, 10.0);
  static std::normal_distribution<double> s1(0, 4.0);
  static std::normal_distribution<double> s2(0, 2.0);
  static std::normal_distribution<double> d0(0, 1.0);
  static std::normal_distribution<double> d1(0, 1.0);
  static std::normal_distribution<double> d2(0, 1.0);

  goal_s[0] += s0(gen);
  goal_s[1] += s1(gen);
  goal_s[2] += s2(gen);

  goal_d[0] += d0(gen);
  goal_d[1] += d1(gen);
  goal_d[2] += d2(gen);

  return pair<State, State>(goal_s, goal_d);
}

Coefficients differentiate(Coefficients coeff) {
  Coefficients new_coeff(coeff.size()-1);
  for (auto i=1; i<coeff.size(); i++) {
    new_coeff[i-1] = (i+1)*coeff[i];
  }
  return new_coeff;
}

double exceeds_max_jerk_cost(Coefficients s_coeff, Coefficients d_coeff, double T, json sensor_fusion) {
  double max_jerk = 10;

  Coefficients s_dot = differentiate(s_coeff);
  Coefficients s_dot_dot = differentiate(s_dot);
  Coefficients jerk_coeff = differentiate(s_dot_dot);

  for (double t=0; t < T; t+=dt) {
    if (fabs(polyeval(jerk_coeff, t)) > max_jerk) {
      return 1.0;
    }
  }

  return 0;
}

double exceeds_max_acceleration_cost(Coefficients s_coeff, Coefficients d_coeff, double T, json sensor_fusion) {
  double max_accel = 10;

  Coefficients s_dot = differentiate(s_coeff);
  Coefficients s_acceleration_coeff = differentiate(s_dot);

  for (double t=0; t < T; t+=dt) {
    if (fabs(polyeval(s_acceleration_coeff, t)) > max_accel) {
      return 1.0;
    }
  }

  return 0;
}

double calculate_cost(Coefficients s_coeff, Coefficients d_coeff, double T, json sensor_fusion) {
  double cost = 0;
  cost += exceeds_max_jerk_cost(s_coeff, d_coeff, T, sensor_fusion);
  cost += exceeds_max_acceleration_cost(s_coeff, d_coeff, T, sensor_fusion);
  return cost;
}

//
// Goal calculate_drive_straight_goal(State start_s, State start_d, double T) {
//   double start_position = start_s[0];
//   double start_speed = start_s[1];
//   double start_acceleration = start_s[2];
//   double max_jerk = 5.0; // meters per second cubed
//   double max_acceleration = 5.0; // meters per second squared
//   double targetSpeed = 20; // meters per second
//
//   double goal_position = 0;
//   double goal_speed = 0;
//   double goal_acceleration = 0;
//
//   // max_jerk
//   // end_acceleration = jerk*time + current_acceleration
//   // time = (end_acceleration-current_acceleration) / jerk;
//   double time_to_max_acceleration = (max_acceleration - start_acceleration) / max_jerk;
//   goal_position += 1.0 / 6.0 * max_jerk * pow(time_to_max_acceleration, 3) + time_to_max_acceleration * start_speed;
//
//   // end_speed = 1 / 2 * jerk*pow(time, 2) + current_speed
//   // time = (current_speed - end_speed) / acceleration
//   double speed_at_max_acceleration = 1.0 / 6.0 * max_jerk*pow(time_to_max_acceleration. 2) + start_speed;
//   double time_to_max_speed = (maxSpeed - speed_at_max_acceleration) / max_acceleration;
//
//   // distance = 1 / 6 * max_jerk*pow(time, 3) + time*current_speed;
//   double minSpeed = - 1.0/2.0 * max_jerk * pow(T, 2) + start_speed;
//   double maxSpeed = 1.0/2.0 * max_jerk * pow(T, 2) + start_speed;
//   cout << minSpeed << " < speed < " << maxSpeed << endl;
//
//   double minDistance = - 1.0/6.0 * max_jerk * pow(T, 3) + T*start_speed;
//   double maxDistance = 1.0/6.0 * max_jerk * pow(T, 3) + T*start_speed;
//   cout << minDistance << " < distance < " << maxDistance << endl;
// }


Path PTG(DriveState state, State start_s, State start_d, json sensor_fusion, int T) {
  cout << "starting PTG" << endl;
  double lane = find_lane(start_d[0]);
  json other_vehicle;
  bool is_following_vehicle = false;
  double max_follow_distance = 40.0;

  // find vehicle to follow
  for (int i=0; i<sensor_fusion.size(); i++) {
    double other_vehicle_s = sensor_fusion[i][5];
    other_vehicle_s -= 4; // measure from back of car

    if (in_same_lane(start_d[0], sensor_fusion[i][6])) {
      // in our lane
      if (other_vehicle_s > start_s[0] && (other_vehicle_s - start_s[0] < max_follow_distance)) {
        other_vehicle = sensor_fusion[i];
        is_following_vehicle = true;
      }
    }
  }

  if (other_vehicle.begin() != other_vehicle.end()) {
    cout << "following vehicle" << endl;
  }

  // double t = T - 4*dt;
  vector<Goal> all_goals;
  // while (t <= T + 4*dt) {
  //   if (state == DriveState::Straight && is_following_vehicle) {
  //
  //   } else if (state == DriveState && !is_following_vehicle) {
  //     Goal goal = calculate_drive_straight_goal();
  //     all_goals.push_back(goal);
  //
  //   }
  //
  //   t += dt;
  // }

  cout << "goals" << endl;

  Coefficients best_s_coeff;
  Coefficients best_d_coeff;
  double lowest_cost = 1e20;

  for (double delta_s=1; delta_s<max_follow_distance; delta_s+= 1) {
    SState goal_s = {start_s[0] + delta_s, 0.0, 0.0};
    DState goal_d = start_d;

    Coefficients s_coeff = JMT(start_s, goal_s, T);
    Coefficients d_coeff = JMT(start_d, goal_d, T);

    double cost = calculate_cost(s_coeff, d_coeff, T, sensor_fusion);

    cout << "start_s: " << start_s[0] << " start_s_generated: " << polyeval(s_coeff, 0) << " end_s: " << goal_s[0] << " cost: " << cost;
    if (cost < lowest_cost) {
      lowest_cost = cost;
      best_s_coeff = s_coeff;
      best_d_coeff = d_coeff;
      cout << " (lowest so far)";
    }

    cout << endl;
    // Coefficients s_coeff = JMT(startS, endS, dt*numPlanningPoints);

  }
  cout << "ending ptg" << endl;
  return Path(best_s_coeff, best_d_coeff, T);
}
