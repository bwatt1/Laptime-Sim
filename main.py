# LIBRARIES
import numpy as np
import matplotlib.pyplot as plt

# PARAMETERS
m = 800 # mass [kg]
rho = 1.2 # air density [kg/m3]
CdA = 1.2 # drag
Crr = 0.015 # rolling resistance coeff
mu = 1.6 # tyre grip
g = 9.81 # law of the universe

P_ice_max = 600e3 # max ICE power by reg [W]
P_ers_max = 120e3 #max ERS power by reg [W]
E_ers_max = 4e6 # max deployable ERS energy by reg [J]

ds = 1.0 # distance step [m]

# TRACK
N = 5000
s = np.arange(0, N*ds, ds)

kappa = np.zeros_like(s)

kappa[1000:1500] = 1/100
kappa[3000:3500] = 1/80

# LATERAL DYNAMICS
v_lat_max = np.sqrt(mu* g /np.maximum(kappa, 1e-6))
v_lat_max[kappa == 0] = 1000 #essentially inf, will be power limited not grip limited here

# PU
def power_unit(v, E_remaining, is_straight):
    P_ice = P_ice_max

    if is_straight and (E_remaining > 0):
        P_ers = P_ers_max
    else: P_ers = 0.0

    P_total = P_ice + P_ers
    return P_total, P_ers


# LONGITUDINAL FORCES
def forces(v):
    F_drag = 0.5 * rho * CdA * v**2
    F_rr = Crr * m * g
    return F_drag, F_rr

# BACKWARDS PASS
v = v_lat_max.copy()

for i in reversed(range(len(s)-1)):
    v_next = v[i+1]
    a_brake = mu * g

    v[i] = min(v[i], np.sqrt(v_next**2 + 2 * a_brake * ds))

# FORWARDS PASS
E_remaining = E_ers_max
for i in range(len(s)-1):
    v_i = v[i]

    is_straight = kappa[i] < 1e-5

    P_total, P_ers = power_unit(v_i, E_remaining, is_straight)

    F_drive = P_total / max(v_i, 1.0)

    F_max = mu * m * g #traction limit
    F_drive = min(F_drive, F_max) #drive force chooses based on power/traction limit

    F_drag, F_rr = forces(v_i)

    a = (F_drive - F_drag - F_rr) / m

    v_next = np.sqrt(max(v_i**2 + 2 * a * ds, 0))

    v[i+1] = min(v[i+1], v_next)

    dt = ds / max(v_i, 1.0)
    E_remaining -= P_ers * dt
    E_remaining = max(E_remaining, 0) # prevent overdeploy

lap_time = np.sum(ds / np.maximum(v, 1.0))
print(lap_time)

# PLOTS
plt.plot(s, v)
plt.xlabel("Distance [m]")
plt.ylabel("Speed [m/s]")
plt.show()

