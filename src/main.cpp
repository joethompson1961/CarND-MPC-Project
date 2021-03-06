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

// Evaluate a polynomial using Horner's method
// Coefficients are ordered lowest order to highest order, e.g. [(x**0), (x**1), ... ,(x**5)]
double polyeval(Eigen::VectorXd coeffs, double x) {
    double result = coeffs[coeffs.size()-1];
    for (size_t i = coeffs.size()-2; i >= 0; i--) {
        result = result*x + coeffs[i];
    }
    return result;
}

//// Evaluate a polynomial.
//double polyeval(Eigen::VectorXd coeffs, double x) {
//  double result = 0.0;
//  for (int i = 0; i < coeffs.size(); i++) {
//    result += coeffs[i] * pow(x, i);
//  }
//  return result;
//}

// Fit a polynomial.
// Adapted from https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int i = 0; i < xvals.size(); i++) {
    for (int j = 0; j < order; j++) {
      A(i, j + 1) = A(i, j) * xvals(i);
    }
  }

  return A.householderQr().solve(yvals);
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

class CmdLineParser{
    public:
        CmdLineParser (int &argc, char **argv){
            for (int i=1; i < argc; ++i)
                this->tokens.push_back(std::string(argv[i]));
        }

        const std::string& getOption(const std::string &cmd) const{
            std::vector<std::string>::const_iterator iterator;
            iterator =  std::find(this->tokens.begin(), this->tokens.end(), cmd);
            if (iterator != this->tokens.end() && ++iterator != this->tokens.end()){
                return *iterator;
            }
            static const std::string empty_string("");
            return empty_string;
        }

        bool exists(const std::string &cmd) const{
            return std::find(this->tokens.begin(), this->tokens.end(), cmd)
                   != this->tokens.end();
        }

    private:
        std::vector<std::string> tokens;
};

void test_poly_fit()
{
    Eigen::VectorXd xvals(25);
    xvals << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24;
    Eigen::VectorXd yvals(25);
    yvals << 1.0, 2.1109, 3.4928, 5.2267, 7.3936, 10.0625, 13.2784, 17.0503, 21.3392, 26.0461, 31.0, 35.9459, 40.5328,
            44.3017, 46.6736, 46.9375, 44.2384, 37.5653, 25.7392, 7.4011, -19.0, -55.2191, -103.2272, -165.2233, -243.6464;
    // Fit a 2nd order polynomial to the trajectory
    auto coeffs = polyfit(xvals, yvals, 5);

    // Result should be:
    //   coeffs:  1,1,0.1,0.01,0.001,-0.0001
    //   for:     f(x) =  1 + x + 0.1*x**2 + 0.01*x**3 + 0.001*x**4 - 0.0001*x**5
    std::cout << "coeffs:" <<coeffs[0] <<"," <<coeffs[1] <<"," <<coeffs[2] <<"," <<coeffs[3] <<"," <<coeffs[4] <<"," <<coeffs[5] <<std::endl;
    exit(0);
}

int main(int argc, char **argv) {

//  test_polyfit();

  CmdLineParser cmdline(argc, argv);
  std::string cmd_str;

  // -n for command line override of N
  int cmd_n = 15;  // default N = 10
  cmd_str = cmdline.getOption("-n");
  if (!cmd_str.empty())
      cmd_n = std::stoi(cmd_str);
  size_t N = cmd_n;

  // -dt for command line override of dt
  double dt = 0.05; // msec
  cmd_str = cmdline.getOption("-dt");
  if (!cmd_str.empty())
    dt = std::stod(cmd_str);

  // Set MPC timestep length, duration, and target velocity
  double ref_v = 125.0;
  cmd_str = cmdline.getOption("-ref_v");
  if (!cmd_str.empty())
    ref_v = std::stod(cmd_str);

  int latency = 100;  // 100 msec actuator latency
  cmd_str = cmdline.getOption("-lat");
  if (!cmd_str.empty())
    latency = std::stoi(cmd_str);

  // -clat for MPC solver latency
  //
  // The latency multiplier adds some additional latency to account for the simulator's processing time plus socket communication delay.
  //
  // total latency = sleep(latency) + MPC solver time
  double clat = 29;
  cmd_str = cmdline.getOption("-clat");
  if (!cmd_str.empty())
    clat = std::stod(cmd_str);

  // -w_cte for command line override of cost function weight for cte
  double w_cte = 70; // cost function weight for cte
  cmd_str = cmdline.getOption("-w_cte");
  if (!cmd_str.empty()){
    w_cte = std::stod(cmd_str);
  }

  double w_epsi   = 1.0;  // psi error
  cmd_str = cmdline.getOption("-w_epsi");
  if (!cmd_str.empty()){
    w_epsi = std::stod(cmd_str);
  }

  double w_ev     = 1.0;   // velocity error
  cmd_str = cmdline.getOption("-w_ev");
  if (!cmd_str.empty()){
    w_ev = std::stod(cmd_str);
  }

  double w_delta  = 1000.0;// steering actuation
  cmd_str = cmdline.getOption("-w_delta");
  if (!cmd_str.empty()){
    w_delta = std::stod(cmd_str);
  }

  double w_a      = 0.1;   // throttle actuation
  cmd_str = cmdline.getOption("-w_a");
  if (!cmd_str.empty()){
    w_a = std::stod(cmd_str);
  }

  double w_sd      = 4500.0;   // rate of steering actuation change
  cmd_str = cmdline.getOption("-w_sd");
  if (!cmd_str.empty()){
    w_sd = std::stod(cmd_str);
  }

  double w_sa      = 50.0;   // rate of throttle actuation change
  cmd_str = cmdline.getOption("-w_sa");
  if (!cmd_str.empty()){
    w_sa = std::stod(cmd_str);
  }

  double w_curve  = 0.12;   // speed around curves
  cmd_str = cmdline.getOption("-w_curve");
  if (!cmd_str.empty()){
    w_curve = std::stod(cmd_str);
  }

  double w_wide  = 1.0;   // wide turns
  cmd_str = cmdline.getOption("-w_wide");
  if (!cmd_str.empty()){
    w_wide = std::stod(cmd_str);
  }

  cout << "N: " << N << endl;     // MPC duration
  cout << "dt: " << dt << endl;   // MPC timestep length
  cout << "ref_v: " << ref_v << endl;  // target velocity
  cout << "latency: " << latency << endl;
  cout << "clat: " << clat << endl;
  cout << "w_cte: " << w_cte << endl;
  cout << "w_epsi: " << w_epsi << endl;
  cout << "w_ev: " << w_ev << endl;
  cout << "w_delta: " << w_delta << endl;
  cout << "w_a: " << w_a << endl;
  cout << "w_sd: " << w_sd << endl;
  cout << "w_sa: " << w_sa << endl;
  cout << "w_curve: " << w_curve << endl;

  // MPC is initialized here!
  MPC mpc(N, dt, ref_v, w_cte, w_epsi, w_ev, w_delta, w_a, w_sd, w_sa, w_curve, w_wide);

  auto start_time = std::chrono::system_clock::now();  // start time
  std::chrono::duration<double> proc_time = std::chrono::duration<double>::zero();
  int ct = 0;

  uWS::Hub h;

  h.onMessage([&mpc, &N, &dt, &clat, &ref_v, &latency, &start_time, &proc_time, &ct](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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
          double psi = j[1]["psi"];  // car heading
          double v = j[1]["speed"];
          double steering = j[1]["steering_angle"];
          double throttle = j[1]["throttle"];

#if 0
          cout << "  car state  x: " << x << "  y: " << y << "  psi: " << psi << "  speed: " << v << endl << endl;
          cout << "  steering: " << steering << "  throttle: " << throttle;
#endif

          // For display purposes only, transform the trajectory to vehicle coordinate system
          // prior to compensating for latency. This keeps the display centered on the road
          // instead of causing it to appear shifted around curves.
          Eigen::VectorXd x_display(ptsx.size());
          Eigen::VectorXd y_display(ptsx.size());
          for (i = 0 ; i < ptsx.size() ; i++) {
              transform_map_to_car(x, y, psi, ptsx[i], ptsy[i], x_display[i], y_display[i]);
          }

          // To compensate for actuation latency, predict the vehicle's state to reflect where it will
          // be after one latency period and use that as the current state.
          double dl = (float)(latency + clat)/1000;  // convert msec to seconds
          x += v * cos(psi) * dl;
          y += v * sin(psi) * dl;
          psi -= v * steering * dl / 2.67;
          v = v + throttle * dl;

          // Transform trajectory waypoints to vehicle's coordinate system (make them relative to the car)
          // The trajectory that we're trying to follow (ptsx[],ptsy[]) is the centerline of the road.
          Eigen::VectorXd x_vals(ptsx.size());
          Eigen::VectorXd y_vals(ptsx.size());
          for (i = 0 ; i < ptsx.size() ; i++) {
              transform_map_to_car(x, y, psi, ptsx[i], ptsy[i], x_vals[i], y_vals[i]);
          }

          // Fit a 2nd order polynomial to the trajectory
          auto coeffs = polyfit(x_vals, y_vals, 3);

          // From world map perspective the cross track error is calculated by evaluating the trajectory
          // polynomial at x, i.e. f(x), and subtracting y, i.e. from where the vehicle is in global frame.
          // The trajectory has been converted to car coordinates (above) so the polynomial
          // is now relative to car and the initial x is 0.  Since the map follows the center of the
          // lane, i.e. y=0 at center, the cte calculation is simplified to be just the evaluation of the
          // polynomial at x = 0.
          double cte = coeffs[0];

          // The orientation error, epsi, is psi - atan(f'(x)) where f'(x) is the slope of the trajectory at position x.
          // Since the polynomial is relative to car then psi = 0 and additionally x = 0 so the calculation of derivative
          // f'(x) is simplified to just coeffs[1] (The coefficients are ordered lowest order to highest order).
          double epsi = -atan(coeffs[1]);  // epsi the difference between car's current heading and the desired heading from trajectory.

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

          auto time1 = std::chrono::system_clock::now();

          vector<double> a1 = mpc.Solve(state, coeffs);

          auto time2 = std::chrono::system_clock::now();
          proc_time += (time2 - time1);

          // Measure average processing time for mpc.Solve()
          ct += 1;
          if (ct == 30)
          {
              ct = 0;
              auto avg = proc_time / 30.0;
              proc_time = std::chrono::duration<double>::zero();
//              std::cout <<"avg:" <<avg.count() <<std::endl;
          }

#if 0
          cout << "a1: " << endl << endl;
          for(i=0; i<a1.size(); ++i)
            std::cout << a1[i] << ' ';
          cout << endl << endl;
#endif

          double steer_value = a1[N*6];
          double throttle_value = a1[N*7-1];

//          cout << "throttle: " << throttle_value << "   steer: " << steer_value << endl;

          /*
           * Build json message for simulator
           */
          json msgJson;
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          // Display the MPC predicted trajectory.
          //
          // Add (x,y) some of the points of the predicted trajectory to list here to have
          // them displayed by the simulator. The points are in reference to the vehicle's
          // coordinate system and are displayed in the simulator connected by a green line.
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;
          for (i = 1 ; i < N ; i++) {
              mpc_x_vals.push_back(a1[0 + i]);
              mpc_y_vals.push_back(a1[N + i]);
          }
          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          // Display the waypoints/reference line.
          //
          // Add (x,y) points of the waypoints to list here to have them displayed by
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
          // The car should be able to drive around the track with 100ms latency.
          this_thread::sleep_for(chrono::milliseconds(latency));

          /*
           * Send json message to simulator
           */
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
  // program doesn't compile :-(
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
