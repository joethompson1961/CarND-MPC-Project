#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// Evaluate a polynomial using Horner's method
// Coefficients are ordered lowest order to highest order, e.g. [(x**0), (x**1), ... ,(x**5)]
AD<double> ADpolyeval(CPPAD_TESTVECTOR(AD<double>) coeffs, AD<double> x) {
    AD<double> result = coeffs[coeffs.size()-1];
    for (int i = coeffs.size()-2; i >= 0; i--) {
        result = result*x + coeffs[i];
    }
    return result;
}

// Returns coefficients for the derivative of a nth order polynomial.
// Coefficients are ordered lowest order to highest order, e.g. [(x**0), (x**1), ... ,(x**5)]
CPPAD_TESTVECTOR(AD<double>) ADdifferentiate(const CPPAD_TESTVECTOR(AD<double>) &coeffs)
{
    CPPAD_TESTVECTOR(AD<double>) new_c(coeffs.size()-1);
    for (size_t n=1; n<coeffs.size(); n++)
        new_c[n-1] = n * coeffs[n];
    return new_c;
}

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to center-of-gravity that has a similar radius.
const double Lf = 2.67;

// The solver receives the state and actuator vectors within one large vector. These indices
// distinguish where within the large vector each variable vector starts and, implicitly, where each ends.
size_t x_start;
size_t y_start;
size_t psi_start;
size_t v_start;
size_t cte_start;
size_t epsi_start;
size_t delta_start;
size_t a_start;

// This class, FG_eval, is a C++ functor, i.e. it overloads operator(), i.e. it returns
// a function that can then be used by the caller. Note that functors are always copied,
// so good to avoid large code and data that is expensive to copy.
class FG_eval {
public:
    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

    ADvector coeffs;  // Fitted polynomial coefficients of the desired trajectory.
    ADvector d_coeffs;  // First derivative of fitted polynomial
    ADvector dd_coeffs;  // Second derivative of fitted polynomial

    size_t N;
    AD<double> dt;
    AD<double> ref_v;

    // Cost function weights with defaults - controls relative cost of each part
    AD<double> w_cte    = 1.0;  // cross track error
    AD<double> w_epsi   = 1.0;  // psi error
    AD<double> w_ev     = 1.0;  // velocity error
    AD<double> w_delta  = 1.0;  // steering actuation
    AD<double> w_a      = 1.0;  // throttle actuation
    AD<double> w_sd     = 1.0;  // rate of steering change
    AD<double> w_sa     = 1.0;  // rate of throttle change
    AD<double> w_curve  = 1.0;  // speed around curves
    AD<double> w_wide   = 1.0;  // wide turns

    // Constructor
    FG_eval(Eigen::VectorXd coeffs, size_t N, double dt, double ref_v, double w_cte, double w_epsi,
            double w_ev, double w_delta, double w_a, double w_sd, double w_sa, double w_curve, double w_wide) {
      this->coeffs.resize(coeffs.size());
      for (int i=0; i<coeffs.size(); i++)
          this->coeffs[i] = coeffs[i];
      this->d_coeffs = ADdifferentiate(this->coeffs);
      this->dd_coeffs = ADdifferentiate(this->d_coeffs);
      this->N = N;
      this->dt = dt;
      this->ref_v = ref_v;

      // cost function weights
      this->w_cte = w_cte;
      this->w_epsi = w_epsi;
      this->w_ev = w_ev;
      this->w_delta = w_delta;
      this->w_a = w_a;
      this->w_sd = w_sd;
      this->w_sa = w_sa;
      this->w_curve = w_curve;
      this->w_wide = w_wide;
    }

    // This functor is used to calculate:
    //   1) What is the cost associated with each set of inputs from the solver.
    //   2) How well the inputs fit within the constraints of the vehicle's motion model.
    // In the general case of optimization solvers this is often referred to as the objective
    // function. It could be used to evaluate any number of additional constraints.  For example
    // it could evaluate the input states against the states of other obstacles looking for
    // collisions, but this isn't applicable to this MPC feedback controller which is used
    // simply to go around the course as fast as possible.
    //
    // In this usage of the solver the input vars represents a sequence of N future vehicle states.
    //
    // The idea here is that the solver is trying to find a set of inputs (i.e. vars - the things
    // that are controlled by the optimizer, in this case the vehicle's states and actuators) that
    // produce the lowest cost while not violating the constraints, i.e. where the contents of all
    // elements in the fg vector ideally converge to zero, indicating a local minimum.
    //
    // For evaluation of constraints, this functor returns the difference between the states
    // chosen/guessed by the solver, i.e. "vars", and the constrained behavior of the vehicle
    // based on the vehicle's motion model.  It does this for each state transition.  In other
    // words, at state N, the motion model predicts what the next state, N+1, should be.  The
    // difference between the guessed next state and the predicted next state are returned.
    // This provides feedback to the solver about how well the guessed states fit the
    // constraints.  Reducing this deviation at all states is the goal of the solver.  It can
    // require many iterations for the solver to converge on a solution.
    //
    // The first and second derivatives of this functor provide feedback to the solver about
    // whether the solution is converging to a local minimum and how fast it's converging.
    // CppAD is used by the solver to generate first and second derivatives. CppAD (Algorithmic
    // Differentiation, aka Automatic Differentiation) is for the step by step conversion from
    // an algorithm that computes function values to an algorithm that computes derivative values.
    // Given a C++ algorithm that computes function values, CppAD generates an algorithm that
    // computes its derivative values.
    void operator()(ADvector& fg, const ADvector& vars) {
        // `vars` is a vector of input variable values (state & actuators), the things the solver
        // controls.
        //
        // `fg` a vector of the cost and constraints to be filled by this functor.
        //
        // The accumulated total cost is stored in the first element of `fg`, i.e. fg[0].
        size_t t;
        AD<double> cost;
        fg[0] = 0;

        // Here we define the costs related to the reference states
        for (t = 0; t < N ; t++) {
            // The radius of curvature at any point x for the curve y = f(x) is given by:
            //    r = (1 + (dy/dx)**2)**(3/2) / (d2y/dx2)
            // The lateral acceleration at any point on the curve is given by:
            //    la = v**2/r
            AD<double> x = vars[x_start + t];
            AD<double> dy_dx = ADpolyeval(this->d_coeffs, x);
            AD<double> d2y_d2x = ADpolyeval(this->dd_coeffs, x);
            AD<double> r = CppAD::pow(1.0 + dy_dx*dy_dx, 1.5) / d2y_d2x;
            AD<double> v = vars[v_start + t];
            AD<double> la = (v*v)/r;

            // Minimize orientation errors
            cost = w_epsi * CppAD::pow(vars[epsi_start + t], 2);
//            cout << "cost - epsi: " << cost << endl;
            fg[0] += cost;

            // Minimize velocity error (the difference between velocity and target velocity)
            cost = w_ev * CppAD::pow(ref_v - vars[v_start + t], 2);
//            cout << "cost - ev: " << cost << endl;
            fg[0] += cost;

            // Minimize cte, but adjust cost based on lateral acceleration to reward taking tight turns.
            AD<double> tgt_cte = vars[cte_start + t] + la * 0.0133; // right turn prefers right side, left turn prefers left side
            cost = w_cte * CppAD::pow(tgt_cte, 2);
            fg[0] += cost;

            // Minimize high lateral accelerations, i.e. going too fast around turns.
            cost = w_curve * CppAD::pow(la, 2);  // always positive
            fg[0] += cost;
        }

        // Minimize the use of actuators.
        for (t = 0; t < N - 1; t++) {
            cost = w_delta * CppAD::pow(vars[delta_start + t], 2);
//            cout << "cost - steering use: " << cost << endl;
            fg[0] += cost;

            cost = w_a * CppAD::pow(vars[a_start + t], 2);
//            cout << "cost - accel use: " << cost << endl;
            fg[0] += cost;
        }

        // Minimize gap between sequential actuations
        for (t = 0; t < N - 2; t++) {
            cost = w_sd * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
            fg[0] += cost;

            cost = w_sa * CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
            fg[0] += cost;
        }

        // Here we setup the model constraints.

        // Adding 1 to each of the starting indices due to cost being located at
        // index 0 of `fg`. This bumps up the position of the other values.
        fg[1 + x_start] = vars[x_start];
        fg[1 + y_start] = vars[y_start];
        fg[1 + psi_start] = vars[psi_start];
        fg[1 + v_start] = vars[v_start];
        fg[1 + cte_start] = vars[cte_start];
        fg[1 + epsi_start] = vars[epsi_start];

        // Set up the rest of the constraints <from the motion model>
        for (t = 1; t < N; t++) {
            AD<double> x0 = vars[x_start + t - 1];
            AD<double> y0 = vars[y_start + t - 1];
            AD<double> psi0 = vars[psi_start + t - 1];
            AD<double> v0 = vars[v_start + t - 1];
            AD<double> cte0 = vars[cte_start + t - 1];
            AD<double> epsi0 = vars[epsi_start + t - 1];
            AD<double> delta0 = vars[delta_start + t - 1];
            AD<double> a0 = vars[a_start + t - 1];
            AD<double> f0 = ADpolyeval(coeffs, x0);
            AD<double> psides0 = CppAD::atan(ADpolyeval(d_coeffs, x0));  // desired psi

            AD<double> x1 = vars[x_start + t];
            AD<double> y1 = vars[y_start + t];
            AD<double> psi1 = vars[psi_start + t];
            AD<double> v1 = vars[v_start + t];
            AD<double> cte1 = vars[cte_start + t];
            AD<double> epsi1 = vars[epsi_start + t];

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
MPC::MPC(size_t N, double dt, double ref_v, double w_cte, double w_epsi,
         double w_ev, double w_delta, double w_a, double w_sd, double w_sa, double w_curve, double w_wide) {
  this->N = N;
  this->dt = dt;
  this->ref_v = ref_v;
  this->w_cte = w_cte;
  this->w_epsi = w_epsi;
  this->w_ev = w_ev;
  this->w_delta = w_delta;
  this->w_a = w_a;
  this->w_sd = w_sd;
  this->w_sa = w_sa;
  this->w_curve = w_curve;
  this->w_wide = w_wide;

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
  FG_eval fg_eval(coeffs, N, dt, ref_v, w_cte, w_epsi, w_ev, w_delta, w_a, w_sd, w_sa, w_curve, w_wide);

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
  options += "Numeric max_cpu_time          0.025\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= (solution.status == CppAD::ipopt::solve_result<Dvector>::success);

  // Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  vector<double> x1;
  for (i = 0; i < n_vars; i++) {
      x1.push_back(solution.x[i]);
  }

  return x1;
}
