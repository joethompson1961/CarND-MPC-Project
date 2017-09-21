#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:
  size_t N;
  double dt;
  double ref_v;
  double w_cte;
  double w_epsi;
  double w_ev;
  double w_delta;
  double w_a;
  double w_sd;
  double w_sa;
  double w_curve;

  MPC(size_t N, double dt, double ref_v, double w_cte, double w_epsi,
      double w_ev, double w_delta, double w_a, double w_sd, double w_sa, double w_curve);

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuations.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
