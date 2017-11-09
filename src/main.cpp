#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include "helpers.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  DriveState state = DriveState::Straight;

  h.onMessage([&state, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;


            tk::spline splineX, splineY, splineDX, splineDY;
            splineX.set_points(map_waypoints_s, map_waypoints_x);
            splineY.set_points(map_waypoints_s, map_waypoints_y);
            splineDX.set_points(map_waypoints_s, map_waypoints_dx);
            splineDY.set_points(map_waypoints_s, map_waypoints_dy);

            // for (int i=0; i<sensor_fusion.size(); i++) {
            //   cout << sensor_fusion[i][0] << " lane " << sensor_fusion[i][6] << endl;
            // }

            cout << "previous path length: " << previous_path_x.size() << endl;

            for (int i=0; i<min((int)previous_path_x.size(), 6); i++) {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
              cout << " x: " << previous_path_x[i] << " y " << previous_path_y[i] << endl;

            }

            int lane_width = 4;
            int lane_1_multiple = lane_width/2;
            int lane_2_multiple = lane_width/2 + lane_width;
            int lane_3_multiple = lane_width/2 + 2*lane_width;

            double numPlanningPoints = 60-next_x_vals.size(); // 0.6 s lookahead
            double dt = 0.02; // 20 ms
            double max_jerk = 10.0;
            // max_jerk
            // acceleration = max_jerk*time
            // speed = 0.5 * max_jerk*pow(time, 2) + current_speed
            // distance = 1 / 6 * max_jerk*pow(time, 3) + time*current_speed;
            double minSpeed = - 1.0/2.0 * max_jerk * pow(numPlanningPoints*dt, 2) + car_speed;
            double maxSpeed = 1.0/2.0 * max_jerk * pow(numPlanningPoints*dt, 2) + car_speed;
            cout << minSpeed << " < speed < " << maxSpeed << endl;

            double minDistance = - 1.0/6.0 * max_jerk * pow(numPlanningPoints*dt, 3) + numPlanningPoints*dt*car_speed;
            double maxDistance = 1.0/6.0 * max_jerk * pow(numPlanningPoints*dt, 3) + numPlanningPoints*dt*car_speed;
            cout << minDistance << " < distance < " << maxDistance << endl;
            double targetSpeed = 20; // meters per second

            double start_dx, start_dy;
            double start_s;
            double start_x = car_x;
            double start_y = car_y;
            double start_x_dot = cos(car_yaw) * car_speed;
            double start_y_dot = sin(car_yaw) * car_speed;
            double start_s_dot;
            if (next_x_vals.size() > 1) {
              start_s = car_s + targetSpeed * next_x_vals.size() * dt;
              start_x = next_x_vals[next_x_vals.size()-1];
              start_y = next_y_vals[next_y_vals.size()-1];
              start_x_dot = (next_x_vals[next_x_vals.size()-1] - next_x_vals[next_x_vals.size()-2])/dt;
              start_y_dot = (next_y_vals[next_y_vals.size()-1] - next_y_vals[next_y_vals.size()-2])/dt;
              start_s_dot = sqrt(start_x_dot*start_x_dot+start_y_dot*start_y_dot);
            } else {
              start_s = car_s;
              start_s_dot = car_speed;
              start_dx = splineDX(car_s);
              start_dy = splineDY(car_s);
            }

            // set goal s
            double goal_s = start_s + max(minDistance, min(maxDistance, targetSpeed * numPlanningPoints * dt));
            for (int i=0; i<sensor_fusion.size(); i++) {
              double other_vehicle_s = sensor_fusion[i][5];
              other_vehicle_s -= 4; // measure from back of car
              double other_vehicle_d = sensor_fusion[i][6];
              if (fabs(other_vehicle_d - lane_2_multiple) <= 2) {
                // in our lane
                if (other_vehicle_s > car_s && other_vehicle_s < goal_s) {
                  goal_s = other_vehicle_s;
                }
              }
            }

            double goal_s_dot = max(minSpeed, min(maxSpeed, targetSpeed));
            double goal_s_dot_dot = 0;
            double goal_dx = splineDX(goal_s);
            double goal_dy = splineDY(goal_s);

            double goal_x = splineX(goal_s) + lane_2_multiple*goal_dx;
            double goal_y = splineY(goal_s) + lane_2_multiple*goal_dy;

            double goal_x_dot = (goal_dx*cos(0.5*M_PI)-goal_dy*sin(0.5*M_PI)) * goal_s_dot;
            double goal_y_dot = (goal_dx*sin(0.5*M_PI)+goal_dy*cos(0.5*M_PI)) * goal_s_dot;

            cout << "start_s: " << start_s << " goal_s: " << goal_s << endl;
            cout << "start_s_dot: " << start_s_dot << " goal_s_dot: " << goal_s_dot << endl;
            cout << "start_x: " << start_x << " goal_x: " << goal_x << endl;
            cout << "start_y: " << start_y << " goal_y: " << goal_y << endl;

            vector<double> startS = { start_s, start_s_dot, 0 };
            vector<double> endS = { goal_s, goal_s_dot, 0 };

            vector<double> startD = { car_d, 0, 0 };

            vector<double> startX = { start_x, start_x_dot, 0 };
            vector<double> endX = { goal_x, goal_x_dot, 0 };

            vector<double> startY = { start_y, start_y_dot, 0 };
            vector<double> endY = { goal_y, goal_y_dot, 0 };

            cout << start_x_dot << " " << start_y_dot << " goal_x_dot " << goal_x_dot << " " << goal_y_dot << endl;

            // Coefficients s_coeff = JMT(startS, endS, dt*numPlanningPoints);
            // Coefficients x_coeff = JMT(startX, endX, dt*numPlanningPoints);
            // Coefficients y_coeff = JMT(startY, endY, dt*numPlanningPoints);

            // cout << "xcoeff" << endl << x_coeff << endl << endl;
            // cout << "ycoeff" << endl << y_coeff << endl << endl;

            // double mid_s = polyeval(s_coeff, 0.5*numPlanningPoints*dt);
            double mid_s = start_s + 0.5*(goal_s-start_s);
            double mid_x = splineX(mid_s) + lane_2_multiple*goal_dx;
            double mid_y = splineY(mid_s) + lane_2_multiple*goal_dy;

            Path path = PTG(state, startS, startD, sensor_fusion, numPlanningPoints*dt);

            VectorXd t_vals(5);
            VectorXd x_vals(3);
            VectorXd y_vals(3);
            VectorXd s_vals(5);
            t_vals << -dt, 0, 0.5*dt*numPlanningPoints, dt*numPlanningPoints, dt+dt*numPlanningPoints;
            x_vals << start_x, mid_x, goal_x;
            y_vals << start_y, mid_y, goal_y;
            s_vals << start_s-car_speed*dt, start_s, mid_s, goal_s, goal_s;
            // Coefficients x_coeff = polyfit(t_vals, x_vals, 1);
            // Coefficients y_coeff = polyfit(t_vals, y_vals, 1);
            Coefficients s_coeff = polyfit(t_vals, s_vals, 2);

            tk::spline splinePartialX, splinePartialY;
            splinePartialX.set_points({start_s, mid_s, goal_s}, {start_x, mid_x, goal_x});
            splinePartialY.set_points({start_s, mid_s, goal_s}, {start_y, mid_y, goal_y});

            double current_speed = car_speed;
            double current_s = start_s;
            double acceleration = 0.0;
            double max_acceleration = 10.0;

            Coefficients path_s_coeff, path_d_coeff;
            double path_T;
            tie(path_s_coeff, path_d_coeff, path_T) = path;

            Coefficients path_s_dot_coeff = differentiate(path_s_coeff);
            Coefficients path_s_dot_dot_coeff = differentiate(path_s_dot_coeff);
            Coefficients path_s_jerk_coeff = differentiate(path_s_dot_dot_coeff);

            for (double t=0; t<path_T; t+= dt) {

                double s = polyeval(path_s_coeff, t);
                double d = polyeval(path_d_coeff, t);
                // double x = splineX(s) + lane_2_multiple*splineDX(s);
                // double y = splineY(s) + lane_2_multiple*splineDY(s);

                double x = splinePartialX(s);
                double y = splinePartialY(s);

                next_x_vals.push_back(x);
                next_y_vals.push_back(y);
                cout << "s " << s << " x: " << x << " y " << y << " speed " << polyeval(path_s_dot_coeff, t) << " a: " << polyeval(path_s_dot_dot_coeff, t) << " jerk: " << polyeval(path_s_jerk_coeff, t) << endl;
            }
            //
            // for (int i=1; i<numPlanningPoints; i++) {
            //   // acceleration = min(max_acceleration, acceleration+dt*max_jerk);
            //   // if (current_speed < targetSpeed) {
            //   //   current_speed += acceleration * dt;
            //   // } else {
            //   //   current_speed -= acceleration * dt;
            //   // }
            //   // current_s += current_speed * dt;
            //
            //   // double x = polyeval(x_coeff, dt*i);
            //   // double y = polyeval(y_coeff, dt*i);
            //   double s = polyeval(s_coeff, dt*i);
            //   // double x = splineX(s) + lane_2_multiple*splineDX(s);
            //   // double y = splineY(s) + lane_2_multiple*splineDY(s);
            //
            //   double x = splinePartialX(s);
            //   double y = splinePartialY(s);
            //
            //   next_x_vals.push_back(x);
            //   next_y_vals.push_back(y);
            //   // cout << "s " << s << " x: " << x << " y " << y << " current_speed " << current_speed << " a: " << acceleration*dt << endl;
            //   cout << s << endl;
            // }
            //
            // cout << endl;
            //
            // int interpolate_count = 100;
            // int lookahead = 1000;
            // int i_start = car_s * interpolate_count + 1;
            // // cout << i_start << " " << i_start+lookahead << endl;
            //
            // // for (int i=i_start; i<i_start+lookahead; i++) {
            // //   int wrapped_i = i % ((map_waypoints_s.size()-1)*interpolate_count);
            // //   int s_i = floor(float(wrapped_i) / interpolate_count);
            // //   int partial_i = wrapped_i % interpolate_count;
            // //   double s0 = map_waypoints_s[s_i];
            // //   double s1 = map_waypoints_s[s_i+1];
            // //   double s = s0 + partial_i * (s1-s0)/interpolate_count;
            // //   next_x_vals.push_back(sX(s)+map_waypoints_dx[s_i]*lane_2_multiple);
            // //   next_y_vals.push_back(sY(s)+map_waypoints_dy[s_i]*lane_2_multiple);
            // //   cout << sX(s)+map_waypoints_dx[s_i]*lane_2_multiple << endl;
            // // }


            for (auto i=0; i<next_x_vals.size()-2; i++) {
              auto dx = next_x_vals[i+1] - next_x_vals[i];
              auto dy = next_y_vals[i+1] - next_y_vals[i];

              auto ax = (next_x_vals[i+2] - next_x_vals[i+1]) - dx;
              auto ay = (next_y_vals[i+2] - next_y_vals[i+1]) - dy;


              if (sqrt(ax*ax+ay*ay) > 0.004) {
                cout << "i: " << i << " ax: " << ax << " ay: " << ay << " dx: " << dx << " dy: " << dy << endl;
              }
            }


          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
