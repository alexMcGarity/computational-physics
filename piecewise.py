import math
import matplotlib.pyplot as plt

g = 9.8  # gravitational acceleration [m/s^2]~label`code:piecewise:defs'
rho = 1.2  # density of air [kg/m^3]
m = 0.45  # mass of soccer ball [kg]
r = 0.11  # radius of soccer ball [m]
cd = 0.5  # typical drag coefficient for sphere

# compute quadratic drag constant [kg/m]
c = rho * math.pi * r**2 * cd / 2

x = 0
y = 0  # starting positions [m]
xlist = []
ylist = []  # initialize position lists
vx = vx0 = 30
vy = vy0 = 20  # starting velocities [m/s]
dt = 0.05  # time step [s] (reduce for accuracy)

# repeat loop while object remains airborn~label`code:piecewise:loop'
while y >= 0:
    # keep list of positions as time passes
    xlist.append(x)
    ylist.append(y)
    # compute velocity magnitude
    v = math.sqrt(vx**2 + vy**2)
    # change in position is velocity * change in time
    x += vx * dt
    y += vy * dt
    # change in velocity is acceleration * change in time
    vx += (-c / m * v * vx) * dt
    vy += (-c / m * v * vy - g) * dt

fig = plt.figure()

# plot trajectory~label`code:piecewise:plot'
times = [dt * i for i in range(1 + int(2 * vy0 / g / dt))]
plt.scatter(
    [vx0 * t for t in times],
    [vy0 * t - g * t**2 / 2 for t in times],
    color="0.5",
    label="No drag",
)
plt.scatter(xlist, ylist, color="0.0", label="c = {:8.2E} [kg/m]".format(c))
plt.title("v0x = {} [m/s], v0y = {} [m/s]".format(vx0, vy0))
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.grid()
plt.legend()
plt.show()

fig.savefig("piecewise.pdf")
