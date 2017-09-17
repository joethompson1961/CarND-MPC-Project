#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return (x) * pi() / 180; }
double rad2deg(double x) { return (x) * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order) {
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

// Transform a point as observed by car (relative to car) to grid map perspective.
void transform_car_to_map(double car_x, double car_y, double car_theta,
                          double observation_x, double observation_y, double& transform_x, double& transform_y) {
#if 0 // DEBUG
    cout << "  observation x: " << observation_x << " y: " << observation_y << endl;
#endif

    transform_x = car_x + (observation_x * cos(car_theta)) - (observation_y * sin(car_theta));
    transform_y = car_y + (observation_x * sin(car_theta)) + (observation_y * cos(car_theta));

    #if 0 // DEBUG
    cout << "  transform   x: " << transform_x << " y: " << transform_y << endl << endl;
#endif
}

// Transform a point from grid map perspective to car's perspective (front of car faces positive x)
void transform_map_to_car(double car_x, double car_y, double car_theta,
                         double map_loc_x, double map_loc_y, double& transform_x, double& transform_y) {
#if 0 // DEBUG
    cout << "  map_loc x: " << map_loc_x << " y: " << map_loc_y << endl;
#endif

    transform_x = (map_loc_x - car_x) * cos(car_theta) + (map_loc_y - car_y) * sin(car_theta);
    transform_y = - (map_loc_x - car_x) * sin(car_theta) + (map_loc_y - car_y) * cos(car_theta);

#if 0 // DEBUG
    cout << "  transform   x: " << transform_x << " y: " << transform_y << endl << endl;
#endif
}

int main() {
  uWS::Hub h;

  // Set MPC timestep length, duration, and target velocity
//  size_t N = 80;
//  double dt = 0.025;   // 50 msec time increment
//  size_t N = 40;
//  double dt = 0.050;   // 25 msec time increment
  size_t N = 20;
  double dt = 0.05;
  double ref_v = 72.0;
  int latency = 100;  // 100 msec latency

  // MPC is initialized here!
  MPC mpc(N, dt, ref_v);

  h.onMessage([&mpc, &N, &dt, &ref_v, &latency](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    size_t i;

    string sdata = string(data).substr(0, length);
#if 0
    cout << sdata << endl;
#endif

    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          // current state
          double x = j[1]["x"];
          double y = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          double steering = j[1]["steering_angle"];
          double throttle = j[1]["throttle"];

#if 0
          cout << endl << "  car state  x: " << x << "  y: " << y << "  psi: " << psi << "  speed: " << v << endl << endl;
#endif

          // For display purposes only, transform the trajectory to vehicle coordinate system
          // prior to compensating for latency.  This keeps the display centered on the road
          // instead of causing it to appear shifted around curves.
          Eigen::VectorXd x_display(ptsx.size());
          Eigen::VectorXd y_display(ptsx.size());
          for (i = 0 ; i < ptsx.size() ; i++) {
              transform_map_to_car(x, y, psi, ptsx[i], ptsy[i], x_display[i], y_display[i]);
          }

          // To compensate for latency, project the vehicle's state to reflect where it will
          // be after one latency period and use that as the current state.
          double dl = (float)latency/1000;  // convert msec to seconds
          dl *= 2.0;
          x += v * cos(psi) * dl;
          y += v * sin(psi) * dl;
          psi -= v * steering * dt / 2.67;
          v = v + throttle * dl;

          // Transform trajectory to vehicle coordinate system (make them relative to the car)
          Eigen::VectorXd x_vals(ptsx.size());
          Eigen::VectorXd y_vals(ptsx.size());
          for (i = 0 ; i < ptsx.size() ; i++) {
              transform_map_to_car(x, y, psi, ptsx[i], ptsy[i], x_vals[i], y_vals[i]);
          }

          // Fit a 2nd order polynomial to the trajectory
          auto coeffs = polyfit(x_vals, y_vals, 2);

          // From world map perspective the cross track error is calculated by evaluating the
          // polynomial at x, i.e. f(x), and subtracting y.  But the trajectory has been
          // converted to car coordinates (above) so the polynomial is now relative to car and
          // the initial x & y are 0.  Hence the cte calclulation is simplified to be just
          // the evaluation of the polynomial at x = 0.
          double cte = polyeval(coeffs, 0);

          // The orientation error, epsi, is psi - f'(x). Since the polynomial is relative to car
          // then psi = 0 and additionally x = 0 so the calculation of derivative f'(x) is simplified
          // to just coeffs[1].
          double epsi = -atan(coeffs[1]);

#if 0
          cout << "cte: " << cte << "  epsi: " << epsi << endl;
#endif

          /*
          * Calculate next predicted steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */
          Eigen::VectorXd state(6);
          state << 0, 0, 0, v, cte, epsi;

          vector<double> a1 = mpc.Solve(state, coeffs);

#if 0
          cout << "a1: " << endl << endl;
          for(i=0; i<a1.size(); ++i)
            std::cout << a1[i] << ' ';
          cout << endl << endl;
#endif

          double steer_value = a1[N*6];
          double throttle_value = a1[N*7-1];

//          cout << "throttle: " << throttle_value << "   steer: " << steer_value << endl;

          json msgJson;
          msgJson["steering_angle"] = -steer_value;
          msgJson["throttle"] = throttle_value;

          // Display the MPC predicted trajectory.
          //
          // Add (x,y) some of the points of the predicted trajectory to list here to have
          // them displayed by the the simulator. The points are in reference to the vehicle's
          // coordinate system and are displayed in the simulator connected by a green line.
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;
          for (i = 0 ; i < N ; i++) {
              mpc_x_vals.push_back(a1[0 + i]);
              mpc_y_vals.push_back(a1[N + i]);
          }
          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          // Display the waypoints/reference line.
          //
          // Add (x,y) points of the reference line to list here to have them displayed by
          // the simulator. The points are in reference to the vehicle's coordinate system
          // and are displayed in the simulator connected by a yellow line.
          vector<double> next_x_vals;
          vector<double> next_y_vals;
          for (int i = 0 ; i < x_vals.size() ; i++) {
            next_x_vals.push_back(x_display[i]);
            next_y_vals.push_back(y_display[i]);
          }
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"steer\"," + msgJson.dump() + "]";

          // Latency
          //
          // The purpose is to mimic real driving conditions where the car
          // doesn't actuate the commands instantly.
          //
          // The car should be to drive around the track with 100ms latency.
          this_thread::sleep_for(chrono::milliseconds(latency));

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
