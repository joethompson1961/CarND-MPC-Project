# CarND-Controls-PID
Self-Driving Car Engineer Nanodegree Program

---

## PID Behavior

A PID controller has 3 elements of control: proportional, differential and integral.

These elements are biased using cooefficients which allows the affects of each element to be increased or decreased to achieve an optimal final outcome, i.e. a controller that is stable over the range of inputs.

Fine tuning the PID controller was largely an empirical problem. By making repeated adjustments to the coefficients and observing the resulting behavior I nudged the controller toward a stable, desirable behavior.

Before fine tuning the PID controller for this problem, i.e. controlling steering steering based on cross track error with car moving at high speed, I evaluated each of the control coefficients (Kp, Ki, Kd) in an ordered way, first focussing on the proportional part of the controller, Kp, while setting the other coefficients to zero. I ran the car at various speeds.  Here's a description of the resulting behaviors:

1. **10 MPH**, Kp = 0.1, Ki = 0.0,  Kd = 0.0:
  - Car navigates the entire track but steering swings very wide around tight curves.  Some steering oscilations present but not too severe.  
2. **20 MPH**, Kp = 0.1, Ki = 0.0,  Kd = 0.0:
  - The first hard turn throws the steering into heavy oscillations, car barely manages to stay on the track but makes a successful loop.
3. **30 MPH**, Kp = 0.1, Ki = 0.0,  Kd = 0.0:
  - Steering oscillations begin immediately and slowly increase heading into first curve where the car loses control and leaves the track.  

Next, holding and Kp and speed at those values, I began incrementing the differential portion of the controller, Kd, to dampen the oscillations in the proportional response.   

4. 30 MPH, Kp = 0.1, Ki = 0.0,  **Kd = 0.1**:
  - Steering oscillations begin immediately but are slightly improved, making it around the first curve before crashing. 
5. 30 MPH, Kp = 0.1, Ki = 0.0,  **Kd = 0.5**:
  - This significantly damped the steering oscillations, significantly reducing the over and under shoot of the proportional part of the controller. The car successfully completes a loop but with several very near misses. 

Holding Kp, Kd and speed at those values I then introduce the integral element, Ki, to the controller, adjusting it in steps as follows:

6. 30 MPH, Kp = 0.1, **Ki = 0.1**,  Kd = 0.5:
  - The steering immediately begins oscillating wildly.  
7. 30 MPH, Kp = 0.1, **Ki = 0.01**,  Kd = 0.5:
  - The steering immediately begins oscillating badly, but improved slightly.  
8. 30 MPH, Kp = 0.1, **Ki = 0.001**,  Kd = 0.5:
  - The steering oscillations are greatly reduced but still there is an acculuation of slowly increasing oscillation that finally causes the car to crash.  
9. 30 MPH, Kp = 0.1, **Ki = 0.0001**,  Kd = 0.5:
  - The vehicle makes a complete loop, sufficient to pass the grader, but steering feels sluglish and more fine tuning necessary to reach higher speeds.  

## PID Tuning Method

In the beginning I manually experimented with various coeffient values to find a set that looked pretty good.  I then increased the target speed and found the results to be less desirable.

To better fine tune the controller at higher speeds I found it important to have an objective method for evaluating the results. I implemented a scoring method by measuring 2 values over approximately 1 lap, i.e. for a fixed number of measurements equivalent to about 1 lap depending on the car's speed:
- **accumulated cross-track error** - the sum of all the absolute values of cross track error at each measurement step.
- **max cross-track error** - the maximum absolute value of cross track error at each measurement step.

The goal is to find the set of coefficients that yield the lowest overall score for both values.

Using this scoring method I increased the speed to 45 MPH and manually tried a variety of different coefficient values folowing a procedure similar to twiddle, recording the scores after each lap.

I observed that the resulting score varied some even when re-run with the same coefficient values and speed; this is possibly due to the simulator having some random element in it's behavioral model. This might have made it difficult to get a good outcome with an automated twiddle but I never tried automating it.

After the best score was achieved I was able to increase the target speed to 55MPH and still navigate the course reliably for many laps. 65MPH was the highest speed tested; it completed a lap but then crashed.
