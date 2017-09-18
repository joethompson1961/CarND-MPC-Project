#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// Evaluate a polynomial.
AD<double> ADpolyeval(Eigen::VectorXd coeffs, AD<double> x) {
    AD<double> result = 0.0;
    for (int i = 0; i < coeffs.size(); i++) {
        result += coeffs[i] * CppAD::pow(x, i);
    }
    return result;
}

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to center-of-gravity that has a similar radius.
const double Lf = 2.67;

// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.
size_t x_start;
size_t y_start;
size_t psi_start;
size_t v_start;
size_t cte_start;
size_t epsi_start;
size_t delta_start;
size_t a_start;

class FG_eval {
public:
    // Fitted polynomial coefficients
    Eigen::VectorXd coeffs;
    size_t N;
    double dt;
    double ref_v;

    FG_eval(Eigen::VectorXd coeffs, size_t N, double dt, double ref_v) {
      this->coeffs = coeffs;
      this->N = N;
      this->dt = dt;
      this->ref_v = ref_v;
    }

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

    void operator()(ADvector& fg, const ADvector& vars) {
        // cost function weights - controls relative cost of each part
        double weights_cte    = 0.0;  // cross track error
        double weights_epsi   = 0.0;  // psi error
        double weights_ev     = 0.0;  // velocity error
        double weights_delta  = 0.0;  // steering actuation
        double weights_a      = 0.0;  // throttle actuation

        weights_cte    = 15.0;  // cross track error
        weights_epsi   = 10.0;  // psi error
        weights_ev     = 1.0;  // velocity error
        weights_delta  = 2200.0;// steering actuation
        weights_a      = 3.0;  // throttle actuation

        size_t t;

        // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
        // NOTE: You'll probably go back and forth between this function and
        // the Solver function below.
        // The cost is stored is the first element of `fg`.
        // Any additions to the cost should be added to `fg[0]`.
        fg[0] = 0;

        // Reference State Cost
        // Define the cost related to the reference state
        // Minimize cross track and orientation errors
        for (t = 0; t < N ; t++) {
            fg[0] += weights_cte * CppAD::pow(vars[cte_start + t], 2);
            fg[0] += weights_epsi * CppAD::pow(vars[epsi_start + t], 2);
        }

        // Minimize velocity error (the difference between velocity and target velocity)
        for (t = 0; t < N ; t++) {
            fg[0] += weights_ev * CppAD::pow(ref_v - vars[v_start + t], 2);
        }

        // Minimize the use of actuators.
        for (t = 0; t < N - 1; t++) {
            fg[0] += weights_delta * CppAD::pow(vars[delta_start + t], 2);
            fg[0] += weights_a * CppAD::pow(vars[a_start + t], 2);
        }

        // Setup the model constraints.
        //
        // Add 1 to each of the starting indices due to cost being located at
        // index 0 of `fg`. This bumps up the position of all the other values.
        fg[1 + x_start] = vars[x_start];
        fg[1 + y_start] = vars[y_start];
        fg[1 + psi_start] = vars[psi_start];
        fg[1 + v_start] = vars[v_start];
        fg[1 + cte_start] = vars[cte_start];
        fg[1 + epsi_start] = vars[epsi_start];

        // Set up the rest of the constraints <from the motion model>>
        for (t = 1; t < N; t++) {
            AD<double> x1 = vars[x_start + t];
            AD<double> y1 = vars[y_start + t];
            AD<double> psi1 = vars[psi_start + t];
            AD<double> v1 = vars[v_start + t];
            AD<double> cte1 = vars[cte_start + t];
            AD<double> epsi1 = vars[epsi_start + t];

            AD<double> x0 = vars[x_start + t - 1];
            AD<double> y0 = vars[y_start + t - 1];
            AD<double> psi0 = vars[psi_start + t - 1];
            AD<double> v0 = vars[v_start + t - 1];
            AD<double> cte0 = vars[cte_start + t - 1];
            AD<double> epsi0 = vars[epsi_start + t - 1];

            AD<double> delta0 = vars[delta_start + t - 1];
            AD<double> a0 = vars[a_start + t - 1];

            AD<double> f0 = ADpolyeval(coeffs, x0);
            AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0);  // desired psi

            // The idea here is to constrain this value to be 0.
            //
            // The use of `AD<double>` and use of `CppAD` is so CppAD can
            // compute derivatives and pass these to the solver.
            fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
            fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
            fg[1 + psi_start + t] = psi1 - (psi0 - v0 * delta0 * dt / Lf);
            fg[1 + v_start + t] = v1 - (v0 + (a0 * dt));
            fg[1 + cte_start + t] = cte1 - (f0 - y0 + v0 * CppAD::sin(epsi0) * dt);
            fg[1 + epsi_start + t] = epsi1 - (psi0 - psides0 + v0 * delta0 * dt / Lf);
        }
    }
};

//
// MPC class definition implementation.
//
MPC::MPC(size_t N, double dt, double ref_v) {
  this->N = N;
  this->dt = dt;
  this->ref_v = ref_v;

  x_start = 0;
  y_start = x_start + N;
  psi_start = y_start + N;
  v_start = psi_start + N;
  cte_start = v_start + N;
  epsi_start = cte_start + N;
  delta_start = epsi_start + N;
  a_start = delta_start + N - 1;
}

MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;

  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  // Set the number of model variables for 6 states and 2 actuators:
  //
  // Note that (N timesteps) ==> (N - 1 actuations)
  size_t n_vars = N * 6 + (N - 1) * 2;

  // Set the number of constraints
  size_t n_constraints = N * 6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (i = 0; i < n_vars; i++) {
    vars[i] = 0.0;
  }

  // Set the initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  // Lower and upper limits for variables.
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // Set all non-actuators upper and lower limits
  // to the max negative and positive values.
  for (i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  for (i = delta_start; i < a_start; i++) {
//    vars_lowerbound[i] = -0.236332;
//    vars_upperbound[i] = 0.236332;
    vars_lowerbound[i] = -0.436332; // -25 degrees
    vars_upperbound[i] = 0.436332;  // 25 degrees
//    vars_lowerbound[i] = -0.52;     // -30 degrees
//    vars_upperbound[i] = 0.52;      // 30 degrees
  }

  // Acceleration/deceleration upper and lower limits.
  for (i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs, N, dt, ref_v);

  // options for IPOPT solver
  std::string options;

  // Uncomment this to print information
  options += "Integer print_level  0\n";
  // Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  vector<double> x1;
  for (i = 0; i < n_vars; i++) {
      x1.push_back(solution.x[i]);
  }

  return x1;
}
