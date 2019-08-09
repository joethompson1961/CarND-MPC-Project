# CarND-MPC-Project
Self-Driving Car Engineer Nanodegree Program

---

## MPC Model

The MPC model used in my implementation includes the following state variables:

**State variables**
- x - the x position of the vehicle in the vehicle's grid map.
- y - the y position of the vehicle in the vehicle's grid map.
- psi - the vehicles orientation relative to the vehicles coordinate system, i.e. front is 0 degrees, increasing counter clockwise. Specified in radians.
- velocity - the vehicle's velocity in the direction of psi. Specified in miles per hour.
- cross-track error
- psi error

**Actuations**
My MPC model predicts the values for two vehicle actuators:
- steering - in range of <-0.43 : 0.43> radians (-25 : 5 degrees).
- throttle - in range of <-1.0 : 1.0>.

**Constraints**
My MPC controller implements a cost function to constrain the vehicles behavior.  My cost function tries to minimize the following:
- cross track error
- psi error
- error in velocity vs. target speed
- changes in steering actuator - used to dampen oscillations.
- changes in throttle actuator - used to dampen oscillations.

My MPC controller also implements a kinematic constraint model for the vehicle's state prediction.  The kinematic constraints are a basic euclidean model. The goal of the kinematic constraints is to minimize changes in the state of the vehicle. 

x1   = x0 + v0 * cos(psi0) * dt					// next x

y1   = y0 + v0 * sin(psi0) * dt					// next y

psi1 = psi0 + v0 * steering * dt / Lf			// next psi

v1   = v0 + (throttle * dt)						// next velocity

cte1 = cte0 + f0 - y0 + v0 * sin(epsi0) * dt    // cross track error

epsi = psi0 - psides0 + v0 * delta0 * dt / Lf	// psi error

## Timestep Length and Elapsed Duration (N & dt)
I tried several combinations of values for timestep length, N, and delta time, dt, to provide a 0.5-2 second horizon for the MPC prediction.

The values for N and dt that I tried ranged from N=<10:50> and dt=<0.02:0.05>. In the end, as I increased speeds and experimented with various cost weights I found that N = 20 and dt = 0.05 produced reliable results with reasonable performance.

## Polynomial Fitting and MPC Preprocessing

Before fitting a polynomial to the waypoints, the trajectory waypoints are converted from the global/map coordinate system to the vehicles's coordinate system. In other words, they are made relative to the vehicle instead of the map.

Doing this transformation simplifies the calculations for cross track error and psi error because initial values for x and y relative to the vehicle are always zero.  Hence the cross track error calculation is simplified to be just the evalution of the polynomial at x = 0, and the psi error calculation is simplfied because the derivitive of f(x) used in the calculation is reduced to just the first coefficient[1].

          // From world map perspective the cross track error is calculated by evaluating the
          // polynomial at x, i.e. f(x), and subtracting y.  But the trajectory has been
          // converted to car coordinates (above) so the polynomial is now relative to car and
          // the initial x & y are 0.  Hence the cte calclulation is simplified to be just
          // the evaluation of the polynomial at x = 0.
          double cte = polyeval(coeffs, 0);

          // The orientation error, epsi, is psi - f'(x). Since the polynomial is relative to car
          // then x = 0 and the calculation of derivative f'(x) is simplified to just coeffs[1].
          double epsi = psi - atan(coeffs[1]);

The waypoint transformation also simplifies the initial state used by MPC, since now the current state x, y and psi are all zero relative to the vehicle.

## Model Predictive Control with Latency

The main loop of the simulator interface (in main.c) implements a 100msec sleep to simulate latency in the response of the actuators as would be typical in a real vehicle.

To compenstate for this latency, the vehicle's state is projected one latency period, i.e. 100 msec, into the future. The compensation is done using the same basic kinetic model calculations used in the MPC solver constraints:

          double dl = (float)latency/1000;  // convert msec to seconds
          x = x + v * cos(psi) * dl + throttle * cos(psi) * dl * dl;
          y = y + v * sin(psi) * dl + throttle * sin(psi) * dl * dl;
          psi = psi + v * steering * dt / 2.67;
          v = v + throttle * dl;
