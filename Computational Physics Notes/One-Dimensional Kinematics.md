+ equations for calculating the movement mechanics for a particle along a single dimension (a line), usually up and down
+ [Wikipedia](https://en.wikipedia.org/wiki/Equations_of_motion) - Go to "Kinematic equations for one particle"
+ In python, we can implement the motion of a particle as a falling object using step-wise integration
	+ First we set the initial conditions for the initial height of the particle, gravitational acceleration, initial velocity, initial time, and our time step
	+ While the height is greater than zero:
		+ update the acceleration (-1 times g with no air resistance)
		+ update the time by adding the time step
		+ update the velocity by adding (a times our time step)
		+ finally update the height by adding (v times our time step)