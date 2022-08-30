import numpy as np
from calculos import *
import environment_engine
import matplotlib.pyplot as plt
import utils

# valores finitos para solução
dx = 0.2
dt = 0.0005
t_max = 10
x_max = 20
y_max = 10
freq = 1

print("started")
array_t, X, Y, array_wave = utils.create_wave(x_max, y_max, t_max, dx, dt)
print("generated wave")
excited_wave = utils.generate_excited_wave(t_max, dt, freq, type="ricker", simu_type = "acoustic")
print("generated excited wave")
environ_params = [
    {
        "type": "base",
        "constant": 0.01
    },
    {
        "type": "borehole",
        "constant": 0.0,
        "x_distance": 1.# 0.1
    }
]
shear_speed, bore_params = environment_engine.create_environment(array_wave, dx, environ_params)#, True)

environ_params = [
    {
        "type": "base",
        "constant": 0.05
    },
    {
        "type": "borehole",
        "constant": 0.5,
        "x_distance": 1.# 0.1
    }
]

pressure_speed = environment_engine.create_environment(array_wave, dx, environ_params)#, True)

lambda_1, mu = utils.generate_mu_lambda(shear_speed, pressure_speed, True)

print("Starting simulation")
dt, array_t, result = solve_elastodynamic_equation(
    initial_condition=array_wave,
    excited_wave=excited_wave,
    steps=[dx, dt],
    mu = mu,
    lambda_1=lambda_1,
    rho=shear_speed,
    max_time=t_max,
    freq=freq,
    borehole_params=bore_params
)
print("Simulation ended")

utils.plot_f_l_frames(result)
utils.plot_response(array_t, result, dt, dx)#, print_freqs=True)
# utils.animate_simulation(array_t, result, 150, file_name="movie.mp4")